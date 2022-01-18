#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "share/grid/user_provided_grids_manager.hpp"

#include "ekat/ekat_assert.hpp"

#include <set>
#include <stdexcept>
#include <string>

namespace scream
{

AtmosphereProcess::AtmosphereProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_comm  (comm)
  , m_params(params)
{}

void AtmosphereProcess::
set_grids_manager (const std::shared_ptr<const GridsManager> grids_manager) {
  m_grids_manager = grids_manager;

  this->set_grids(m_grids_manager);
}

void AtmosphereProcess::initialize (const TimeStamp& t0, const RunType run_type) {
  set_fields_and_groups_pointers();
  m_time_stamp = t0;
  initialize_impl(run_type);

  if (m_params.isParameter("Number of Subcycles")) {
    m_num_subcycles = m_params.get<int>("Number of Subcycles");
  }
  EKAT_REQUIRE_MSG (m_num_subcycles>0,
      "Error! Invalid number of subcycles.\n"
      "  - Atm proc name: " + this->name() + "\n"
      "  - Num subcycles: " + std::to_string(m_num_subcycles) + "\n");

  bool do_out_req_fields  = m_params.isSublist("Output Requied Fields");
  bool do_out_comp_fields = m_params.isSublist("Output Computed Fields");
  if (do_out_req_fields || do_out_comp_fields) {
    // Create helper fields manager, adding all fields/groups used by this process
    std::map<std::string,std::shared_ptr<FieldManager<Real>>> field_mgrs;
    for (const auto& gn : this->get_required_grids()) {
      field_mgrs.emplace(gn,m_grids_manager->get_grid(gn));
      field_mgrs[gn]->registration_begins();
      field_mgrs[gn]->registration_ends();
    }

    std::map<std::string,std::vector<std::string>> f_in, f_out;
    for (const auto& f : get_fields_in()) {
      const auto& gn = f.get_header().get_identifier().get_grid_name();
      field_mgrs.at(gn)->add_field(f);
      f_in[gn].push_back(f.get_header().get_identifier().name());
    }
    for (const auto& f : get_fields_out()) {
      const auto& gn = f.get_header().get_identifier().get_grid_name();
      field_mgrs.at(gn)->add_field(f);
      f_out[gn].push_back(f.get_header().get_identifier().name());
    }
    for (const auto& g : get_groups_in()) {
      if (g.m_bundle) {
        auto f = *g.m_bundle;
        const auto& gn = g.grid_name();
        field_mgrs.at(gn)->add_field(f);
        f_in[gn].push_back(f.get_header().get_identifier().name());
      }
    }
    for (const auto& g : get_groups_out()) {
      if (g.m_bundle) {
        auto f = *g.m_bundle;
        const auto& gn = g.grid_name();
        field_mgrs.at(gn)->add_field(f);
        f_out[gn].push_back(f.get_header().get_identifier().name());
      }
    }

    // Helper lambda, to set defaults and mandatory values for output manager param lists
    auto set_defaults = [&](ekat::ParameterList& pl, bool inputs) {
      pl.set<std::string>("Averaging Type", "Instant");
      auto& out_control = pl.sublist("Output Contol");
      out_control.set("Frequency",1);
      out_control.set("Frequency Units","Steps");
      auto& checkpt_control = pl.sublist("Checkpoint Control");
      checkpt_control.set("Frequency",0);

      // Add fields names only if not already present
      if (not pl.isParameter("Field Names")) {
        const auto& grids = this->get_required_grids();
        auto& default_f_names = inputs ? f_in : f_out;
        for (const auto& gn : grids) {
          if (not pl.sublist("Fields").isSublist(gn)) {
            auto fn_pl = pl.sublist("Fields").sublist(gn);
            fn_pl.set<std::vector<std::string>>("Field Names",default_f_names.at(gn));
          }
        }
      }
    };

    // Create output managers for in and/or out fields
    // Note: even if this is a reestarted run, we do not have the possibility of
    //       retrieving the "case" t0, so we just pass 'false' for is_restarted_run.
    if (do_out_req_fields) {
      m_output_required_fields = std::make_shared<OutputManager>();
      auto& out_req_fields_pl = m_params.sublist("Output Requied Fields");
      set_defaults(out_req_fields_pl,true);
      m_output_required_fields->setup(this->get_comm(),out_req_fields_pl,field_mgrs,m_grids_manager,t0,false,false);
    }
    if (do_out_comp_fields) {
      m_output_computed_fields = std::make_shared<OutputManager>();
      auto& out_comp_fields_pl = m_params.sublist("Output Computed Fields");
      set_defaults(out_comp_fields_pl,false);
      m_output_computed_fields->setup(this->get_comm(),out_comp_fields_pl,field_mgrs,m_grids_manager,t0,false,false);
    }
  }
}

void AtmosphereProcess::run (const int dt) {
  if (m_params.get("Enable Input Field Checks", true)) {
    // Run any check on required fields that has been stored in this AP
    check_required_fields();
  }

  EKAT_REQUIRE_MSG ( (dt % m_num_subcycles)==0,
      "Error! The number of subcycle iterations does not exactly divide the time step.\n"
      "  - Atm proc name: " + this->name() + "\n"
      "  - Num subcycles: " + std::to_string(m_num_subcycles) + "\n"
      "  - Time step    : " + std::to_string(dt) + "\n");

  if (m_output_required_fields) {
    m_output_required_fields->run(m_time_stamp+dt);
  }

  // Let the derived class do the actual run
  auto dt_sub = dt / m_num_subcycles;
  for (int isub=0; isub<m_num_subcycles; ++isub) {
    run_impl(dt_sub);
  }

  if (m_params.get("Enable Output Field Checks", true)) {
    // Run any check on required fields that has been stored in this AP
    check_computed_fields();
  }
  if (m_output_computed_fields) {
    m_output_computed_fields->run(m_time_stamp+dt);
  }

  // Update all output fields time stamps
  m_time_stamp += dt;
  update_time_stamps ();
}

void AtmosphereProcess::finalize (/* what inputs? */) {
  finalize_impl(/* what inputs? */);
}

void AtmosphereProcess::set_required_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_field(f.get_header().get_identifier()),
    "Error! Input field is not required by this atm process.\n"
    "    field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "    atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_fields_in,f)) {
    m_fields_in.emplace_back(f);
  }

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as customer if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as customer to the field
    add_me_as_customer(f);
  }

  set_required_field_impl (f);
}

void AtmosphereProcess::set_computed_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_field(f.get_header().get_identifier()),
    "Error! Input field is not computed by this atm process.\n"
    "   field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_fields_out,f)) {
    m_fields_out.emplace_back(f);
  }

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as provider if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as provider for the field
    add_me_as_provider(f);
  }

  set_computed_field_impl (f);
}

void AtmosphereProcess::set_required_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not require the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_groups_in,group)) {
    m_groups_in.emplace_back(group);
    // AtmosphereProcessGroup is just a "container" of *real* atm processes,
    // so don't add me as customer if I'm an atm proc group.
    if (this->type()!=AtmosphereProcessType::Group) {
      if (group.m_bundle) {
        add_me_as_customer(*group.m_bundle);
      } else {
        for (auto& it : group.m_fields) {
          add_me_as_customer(*it.second);
        }
      }
    }
  }

  set_required_group_impl(group);
}

void AtmosphereProcess::set_computed_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not compute the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_groups_out,group)) {
    m_groups_out.emplace_back(group);
    // AtmosphereProcessGroup is just a "container" of *real* atm processes,
    // so don't add me as provider if I'm an atm proc group.
    if (this->type()!=AtmosphereProcessType::Group) {
      if (group.m_bundle) {
        add_me_as_provider(*group.m_bundle);
      } else {
        for (auto& it : group.m_fields) {
          add_me_as_provider(*it.second);
        }
      }
    }
  }

  set_computed_group_impl(group);
}

void AtmosphereProcess::check_required_fields () const { 
  // Loop over all the field checks on input fields, and execute them
  for (const auto& fp_it : m_property_checks_in) {
    const auto& fid   = fp_it.first;
    const auto& fname = fid.name();
    const auto& gname = fid.get_grid_name();
    const auto& prop_check = *fp_it.second;
    if (has_required_field(fid)) {
      auto& f = get_field_in(fname,gname);
      EKAT_REQUIRE_MSG(prop_check.check(f),
          "Error! Input field property check failed.\n"
          "       Atm proc: " + this->name() + "\n"
          "       Field:    " + fname + "\n"
          "       Check:    " + prop_check.name() + "\n"
          " NOTE: we don't allow repairing input fields.\n");
    } else {
      auto& group = get_group_in(fname,gname);
      if (group.m_bundle) {
        EKAT_REQUIRE_MSG(prop_check.check(*group.m_bundle),
            "Error! Input field property check failed.\n"
            "       Atm proc: " + this->name() + "\n"
            "       Field:    " + group.m_info->m_group_name + "\n"
            "       Check:    " + prop_check.name() + "\n"
            " NOTE: we don't allow repairing input fields.\n");
      } else {
        for (const auto& it : group.m_fields) {
          const auto& f = *it.second;
          EKAT_REQUIRE_MSG(prop_check.check(f),
              "Error! Input field property check failed.\n"
              "       Atm proc: " + this->name() + "\n"
              "       Field:    " + fname + "\n"
              "       Check:    " + prop_check.name() + "\n"
              " NOTE: we don't allow repairing input fields.\n");
        }
      }
    }
  }
}

void AtmosphereProcess::check_computed_fields () {
  // Loop over all the field checks on input fields, and execute them
  for (const auto& fp_it : m_property_checks_out) {
    const auto& fid   = fp_it.first;
    const auto& fname = fid.name();
    const auto& gname = fid.get_grid_name();
    const auto& prop_check = *fp_it.second;
    if (has_computed_field(fid)) {
      auto& f = get_field_out(fname,gname);
      const auto check = prop_check.check(f);
      EKAT_REQUIRE_MSG(check || prop_check.can_repair(),
          "Error! Output field property check failed (and cannot be repaired).\n"
          "       Atm proc: " + this->name() + "\n"
          "       Field:    " + fname + "\n"
          "       Check:    " + prop_check.name() + "\n");
      if (not check) {
        prop_check.repair(f);
      }
    } else {
      auto& group = get_group_out(fname,gname);
      if (group.m_bundle) {
        auto& f = *group.m_bundle;
        const auto check = prop_check.check(f);
        EKAT_REQUIRE_MSG(check || prop_check.can_repair(),
            "Error! Output field property check failed (and cannot be repaired).\n"
            "       Atm proc: " + this->name() + "\n"
            "       Field:    " + group.m_info->m_group_name + "\n"
            "       Check:    " + prop_check.name() + "\n");
        if (not check) {
          prop_check.repair(f);
        }
      } else {
        for (const auto& it : group.m_fields) {
          auto& f = *it.second;
          const auto check = prop_check.check(f);
          EKAT_REQUIRE_MSG(check || prop_check.can_repair(),
              "Error! Output field property check failed (and cannot be repaired).\n"
              "       Atm proc: " + this->name() + "\n"
              "       Field:    " + it.first + "\n"
              "       Check:    " + prop_check.name() + "\n");
          if (not check) {
            prop_check.repair(f);
          }
        }
      }
    }
  }
}

bool AtmosphereProcess::has_required_field (const FieldIdentifier& id) const {
  for (const auto& it : m_required_field_requests) {
    if (it.fid==id) {
      return true;
    }
  }
  return false;
}
bool AtmosphereProcess::has_computed_field (const FieldIdentifier& id) const {
  for (const auto& it : m_computed_field_requests) {
    if (it.fid==id) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_required_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_required_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}
bool AtmosphereProcess::has_computed_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_computed_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}

void AtmosphereProcess::update_time_stamps () {
  const auto& t = timestamp();

  // Update *all* output fields/groups, regardless of whether
  // they were touched at all during this time step.
  // TODO: this might have to be changed
  for (auto& f : m_fields_out) {
    f.get_header().get_tracking().update_time_stamp(t);
  }
  for (auto& g : m_groups_out) {
    if (g.m_bundle) {
      g.m_bundle->get_header().get_tracking().update_time_stamp(t);
    } else {
      for (auto& f : g.m_fields) {
        f.second->get_header().get_tracking().update_time_stamp(t);
      }
    }
  }
}

void AtmosphereProcess::add_me_as_provider (const Field& f) {
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

void AtmosphereProcess::add_me_as_customer (const Field& f) {
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void AtmosphereProcess::add_internal_field (const Field& f) {
  m_internal_fields.push_back(f);
}

const Field& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) const {
  return get_field_in_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) {
  return get_field_in_impl(field_name,grid_name);
}

const Field& AtmosphereProcess::
get_field_in(const std::string& field_name) const {
  return get_field_in_impl(field_name);
}

Field& AtmosphereProcess::
get_field_in(const std::string& field_name) {
  return get_field_in_impl(field_name);
}

const Field& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) const {
  return get_field_out_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) {
  return get_field_out_impl(field_name,grid_name);
}

const Field& AtmosphereProcess::
get_field_out(const std::string& field_name) const {
  return get_field_out_impl (field_name);
}

Field& AtmosphereProcess::
get_field_out(const std::string& field_name) {
  return get_field_out_impl (field_name);
}

const FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) const {
  return get_group_in_impl (group_name,grid_name);
}

FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) {
  return get_group_in_impl (group_name,grid_name);
}

const FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name) const {
  return get_group_in_impl(group_name);
}

FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name) {
  return get_group_in_impl(group_name);
}

const FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) const {
  return get_group_out_impl(group_name,grid_name);
}

FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) {
  return get_group_out_impl(group_name,grid_name);
}

const FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name) const {
  return get_group_out_impl(group_name);
}

FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name) {
  return get_group_out_impl(group_name);
}

const Field& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) const {
  return get_internal_field_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) {
  return get_internal_field_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_internal_field(const std::string& field_name) {
  return get_internal_field_impl(field_name);
}

const Field& AtmosphereProcess::
get_internal_field(const std::string& field_name) const {
  return get_internal_field_impl(field_name);
}

void AtmosphereProcess::set_fields_and_groups_pointers () {
  for (auto& f : m_fields_in) {
    const auto& fid = f.get_header().get_identifier();
    m_fields_in_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
  for (auto& f : m_fields_out) {
    const auto& fid = f.get_header().get_identifier();
    m_fields_out_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
  for (auto& g : m_groups_in) {
    const auto& group_name = g.m_info->m_group_name;
    m_groups_in_pointers[group_name][g.grid_name()] = &g;
  }
  for (auto& g : m_groups_out) {
    const auto& group_name = g.m_info->m_group_name;
    m_groups_out_pointers[group_name][g.grid_name()] = &g;
  }
  for (auto& f : m_internal_fields) {
    const auto& fid = f.get_header().get_identifier();
    m_internal_fields_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
}

void AtmosphereProcess::
alias_field_in (const std::string& field_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_fields_in_pointers[alias_name][grid_name] = m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field for aliasing request.\n"
        "    - field name: " + field_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_field_out (const std::string& field_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  try {
    m_fields_out_pointers[alias_name][grid_name] = m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field for aliasing request.\n"
        "    - field name: " + field_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_group_in (const std::string& group_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_groups_in_pointers[alias_name][grid_name] = m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_group_out (const std::string& group_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  try {
    m_groups_out_pointers[alias_name][grid_name] = m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_in_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_in_impl(const std::string& field_name) const {
  try {
    auto& copies = m_fields_in_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find input field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_out_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_out_impl(const std::string& field_name) const {
  try {
    auto& copies = m_fields_out_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find output field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_in_impl(const std::string& group_name, const std::string& grid_name) const {
  try {
    return *m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_in_impl(const std::string& group_name) const {
  try {
    auto& copies = m_groups_in_pointers.at(group_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find input group providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  group name: " + group_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_out_impl(const std::string& group_name, const std::string& grid_name) const {
  try {
    return *m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_out_impl(const std::string& group_name) const {
  try {
    auto& copies = m_groups_out_pointers.at(group_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find output group providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  group name: " + group_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
}

Field& AtmosphereProcess::
get_internal_field_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_internal_fields_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate internal field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_internal_field_impl(const std::string& field_name) const {
  try {
    auto& copies = m_internal_fields_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find internal field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate internal field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

} // namespace scream
