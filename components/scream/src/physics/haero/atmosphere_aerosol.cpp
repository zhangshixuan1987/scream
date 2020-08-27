#include "physics/haero/atmosphere_aerosol.hpp"
#include "physics/haero/haero_inputs_initializer.hpp"

#include "ekat/ekat_assert.hpp"

namespace scream {

/*
  HAero -- Hi-res AEROsol routines
*/

// =========================================================================================
HiResAerosol::HiResAerosol(const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_haero_comm(comm)
  , m_haero_params(params)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, haero options.
*/

  if (!m_haero_params.isParameter("Grid")) {
    m_haero_params.set("Grid", std::string("SE Physics"));
  }

  m_initializer = create_field_initializer<HaeroInputsInitializer>();
}
// =========================================================================================

void HiResAerosol::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  constexpr int NVL = 72;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  constexpr int QSZ =  35;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  constexpr int NMODES = 4;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  using namespace ekat::units;
  auto mass_mixing_ratio_units = kg/kg;
  mass_mixing_ratio_units.set_string("kg-aerosol/kg-air");
  auto number_mixing_ratio_units = 1/kg;
  number_mixing_ratio_units.set_string("#/kg");
  auto conversion_rate_units = mol/(s * mol);
  conversion_rate_units.set_string("mol-new/(mol-old * s)");

  const auto& grid_name = m_haero_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  const int num_dofs = grid->get_num_local_dofs();

  using namespace ShortFieldTagsNames;
  FieldLayout scalar3d_layout_mid {{COL,VL}, {nc, NVL}};
  FieldLayout scalar3d_layout_int {{COL,VL}, {nc, NVL+1}};
  FieldLayout tracers_layout {{COL, VAR, VL}, {nc, QSZ, NVL}};
  FieldLayout aerosol_layout {{COL, VAR, VL, AEROM}, {nc, QSZ, NVL, NMODES}};

  // inputs
  m_required_fields.emplace("T",    scalar3d_layout_mid, K,  grid_name);
  m_required_fields.emplace("pmid", scalar3d_layout_mid, Pa, grid_name);
  m_required_fields.emplace("zi",   scalar3d_layout_int, m,  grid_name);

  // input-outputs
  m_required_fields.emplace("num_aer", aerosol_layout, number_mixing_ratio_units, grid_name);
  m_computed_fields.emplace("num_aer", aerosol_layout, number_mixing_ratio_units, grid_name);

  m_required_fields.emplace("newnuc_soa", aerosol_layout, conversion_rate_units, grid_name);
  m_computed_fields.emplace("newnuc_soa", aerosol_layout, conversion_rate_units, grid_name);

  m_required_fields.emplace("soa_aer", aerosol_layout, mass_mixing_ratio_units, grid_name);
  m_computed_fields.emplace("soa_aer", aerosol_layout, mass_mixing_ratio_units, grid_name);
}
// =========================================================================================

void HiResAerosol::initialize(const util::Timestamp& t0)
{
  m_current_ts = t0;
}
// =========================================================================================


} // namespace scream
