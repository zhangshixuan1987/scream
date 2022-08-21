#include "scream_scorpio_interface.hpp"
#include "scream_scorpio_interface_details.hpp"

#include "share/scream_types.hpp"
#include "scream_config.h"

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <pio.h>

#ifdef SCREAM_CIME_BUILD
extern "C"
{
int get_io_sys        (const int atm_id);
int get_io_type       (const int atm_id);
int get_io_rearranger (const int atm_id);
int get_io_format     (const int atm_id);
}
#endif

namespace scream {
namespace scorpio {

template<typename T>
std::string type2str ();

template<> std::string type2str<int> () { return "int"; }
template<> std::string type2str<float> () { return "float"; }
template<> std::string type2str<double> () { return "double"; }
template<> std::string type2str<std::string> () { return "char"; }

inline std::string e2str (const FileMode mode) {
  return mode==FileMode::Read ? "Read" : "Write";
}

// ============================================================== //

// ------------------------ IMPLEMENTATION ------------------------- //

void init(const ekat::Comm& comm, const int atm_id) {
  auto& s = IOSession::instance(false);

  s.comm = comm;
#ifdef SCREAM_CIME_BUILD
  // PIO system should alreay be inited. Grab it.
  s.iosysid = get_io_sys(atm_id);
  s.iotype  = get_io_type(atm_id);
  s.rearr   = get_io_rearranger(atm_id);
  s.format  = get_io_format(atm_id);
#else
#ifdef PIO_USE_PNETCDF 
  s.iotype = PIO_IOTYPE_PNETCDF;
#else
  s.iotype = PIO_IOTYPE_NETCDF4P;
#endif
  s.rearr  = PIO_REARR_SUBSET;
  s.format = PIO_64BIT_DATA;

  const int stride = 1;
  const int base   = 0;
  int err = PIOc_Init_Intracomm(comm.mpi_comm(),comm.size(),stride,base,s.rearr,&s.iosysid);

  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Could not initialize pio subsystem.\n"
      " - atm rank: " + std::to_string(comm.rank()) + "\n"
      " - pio err : " + std::to_string(err) + "\n");

  s.must_finalize = true;
#endif
}

void finalize() {
  auto& s = IOSession::instance();

  if (s.must_finalize) {
    int err = PIOc_finalize(s.iosysid);

    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Could not finalize pio subsystem.\n"
        " - atm rank: " + std::to_string(s.comm.rank()) + "\n"
        " - pio err : " + std::to_string(err) + "\n");

    s.iosysid = -1;
  }
}

void open_file (const std::string& filename, const FileMode mode) {
  auto& s = IOSession::instance();

  auto& f = s.files[filename];

  if (f!=nullptr) {
    // The file was already open. Make sure it's a read operation
    EKAT_REQUIRE_MSG (mode==f->mode && mode==FileMode::Read,
        "Error! Multiple access to IO file only supported in read mode.\n"
        " - file name: " + filename + "\n"
        " - curr mode: " + e2str(f->mode) + "\n" 
        " - new  mode: " + e2str(mode) + "\n");
  } else {
    f = std::make_shared<IOFile>();
    f->name = filename;
    f->mode = mode;

    int ierr;
    if (mode==FileMode::Read) {
      ierr = PIOc_openfile(s.iosysid,&f->ncid,&s.iotype,filename.c_str(),PIO_NOWRITE);
      f->load_all_entities_from_file ();
      f->set_vars_decomps(s.decomps);
      f->enddef = true;
    } else {
      ierr = PIOc_createfile(s.iosysid,&f->ncid,&s.iotype,filename.c_str(),PIO_WRITE | PIO_CLOBBER);
      register_dimension(filename,"time",PIO_UNLIMITED);
    }
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not open pio file.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    // Immediately register
  }

  ++f->customers;
}
/* ----------------------------------------------------------------- */
void release_file (const std::string& filename) {
  auto& s = IOSession::instance();

  EKAT_REQUIRE_MSG (s.files.find(filename)!=s.files.end(),
      "Error! Could not close file, since it was not found.\n"
      " - file name: " + filename + "\n");

  auto& f = s.get_file(filename);
  --f.customers;

  if (f.customers==0) {
    // Remove vars as customers of their decomp
    for (auto it : f.vars) {
      auto& v = *it.second;
      if (v.decomp!=nullptr) {
        --v.decomp->customers;
      }
    }

    int ierr = PIOc_sync(f.ncid);
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not sync file.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    ierr = PIOc_closefile(f.ncid);
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not close file.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    // File is no longer needed. Clean it up
    s.files.erase(filename);
  }
}
/* ------------------------------------------------- */
void register_dimension (const std::string& filename,
                         const std::string& dim_name,
                         const int length)
{
  auto& s = IOSession::instance();

  int ierr;
  auto& f = s.get_file(filename);
  auto& d = f.dims[dim_name];
  if (d==nullptr) {
    EKAT_REQUIRE_MSG (f.mode==FileMode::Write,
        "Error! Registering new dimension in a read-only file!\n");
    d = std::make_shared<IOEntity>();
    d->name = dim_name;
    ierr = PIOc_def_dim (f.ncid,dim_name.c_str(),length,&d->ncid);
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not define dimension.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - dim name : " + dim_name + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");
  } else {
    // Check that the dimension specs are the same
    const auto name = detail::get_dim_name(f.ncid,d->ncid);
    const auto len  = detail::get_dim_len(f.ncid,d->ncid);

    EKAT_REQUIRE_MSG (name==dim_name,
        "Error! Something is amiss with dimensions storage.\n"
        "       Please, contact developers.\n");
    EKAT_REQUIRE_MSG (len==length,
        "Error! Dimension already registered with different length.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - dim name : " + dim_name + "\n"
        " - dim len  : " + std::to_string(len) + "\n"
        " - len in   : " + std::to_string(length) + "\n");
  }
}
/* ----------------------------------------------------------------- */
void register_variable (const std::string& filename,
                        const std::string& varname,
                        const std::string& longname,
                        const std::string& units,
                        const std::vector<std::string>& dim_names,
                        const std::string& dtype)
{
  auto& s = IOSession::instance();

  int ierr;
  auto& f = s.get_file(filename);
  auto& v = f.vars[varname];
  if (v==nullptr) {
    v = std::make_shared<IOVar>();
    v->name = varname;
    EKAT_REQUIRE_MSG (f.mode==FileMode::Write,
        "Error! Registering new variable in a read-only file!\n");

    std::vector<int> dimids;
    if (dim_names[0]!="time") {
      dimids.push_back(f.get_dim("time").ncid);
    }
    for (const auto& dn : dim_names) {
      dimids.push_back(f.get_dim(dn).ncid);
    }
    ierr = PIOc_def_var(f.ncid, varname.c_str(), detail::str2nctype(dtype),
                        dimids.size(), dimids.data(), &v->ncid);

    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not define variable.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - var name : " + varname + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    // Add units/longname as attribute
    if (units!="") {
      set_attribute(filename,varname,"units",units);
    }
    if (longname!="") {
      set_attribute(filename,varname,"long_name",longname);
    }
  } else {
    // Check that var specs are the same
    EKAT_REQUIRE_MSG (detail::get_var_name(f.ncid,v->ncid)==varname,
        "Error! Something is amiss with dimensions storage.\n"
        "       Please, contact developers.\n");
    const auto nctype = detail::get_var_dtype(f.ncid,v->ncid);
    EKAT_REQUIRE_MSG (nctype==detail::str2nctype(dtype),
        "Error! Variable already defined with different data type.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - var name : " + varname + "\n"
        " - var type : " + detail::nctype2str(nctype) + "\n"
        " - type in  : " + dtype + "\n");
    EKAT_REQUIRE_MSG ("time"+ekat::join(dim_names,"")==ekat::join(detail::get_var_dim_names(f.ncid,v->ncid),""),
        "Error! Variable already defined with different dimensions.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - var name : " + varname + "\n"
        " - var dims : (" + ekat::join(detail::get_var_dim_names(f.ncid,v->ncid),", ") + ")\n"
        " - dims in  : (" + ekat::join(dim_names,", ") + ")\n");
  }
}
/* ----------------------------------------------------------------- */
void register_decomp (const std::string& dtype,
                      const std::vector<std::string>& dims,
                      const std::vector<int>& gdimlen,
                      const std::vector<offset_t>& my_offsets)
{
  auto& s = IOSession::instance();

  EKAT_REQUIRE_MSG (dims.size()==gdimlen.size(),
      "[register_decomp] Error! Input arrays have different lengths.\n"
      " - dims size   : " + std::to_string(dims.size()) + "\n"
      " - gdimlen size: " + std::to_string(gdimlen.size()) + "\n");

  int ierr;
  // convert dtype->nctype->dtype, to get rid of 'real'
  // and use the actual dtype (float or double)
  const auto decomp_name = detail::nctype2str(detail::str2nctype(dtype)) + "."
                         + ekat::join(dims,"-") + "."
                         + ekat::join(gdimlen,"-");

  auto& decomp = s.decomps[decomp_name];
  if (decomp==nullptr) {
    // Create decomp
    decomp = std::make_shared<IODecomp>();
    decomp->name = decomp_name;
    decomp->size = my_offsets.size();

    auto offsets = reinterpret_cast<const PIO_Offset*>(my_offsets.data());
    ierr = PIOc_init_decomp (s.iosysid, detail::str2nctype(dtype),
                             dims.size(), gdimlen.data(),
                             my_offsets.size(), offsets,
                             &decomp->ncid, s.rearr,
                             nullptr, nullptr);

    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not add decomposition.\n"
        " - atm rank: " + std::to_string(s.comm.rank()) + "\n"
        " - decomp  : " + decomp_name + "\n"
        " - pio err : " + std::to_string(ierr) + "\n");
  } else {
    // Check same decomp specs (not much we can check; only size)
    EKAT_REQUIRE_MSG (decomp->size==static_cast<PIO_Offset>(my_offsets.size()),
        "Error! Decomp already registered with different size.\n"
        " - decomp  : " + decomp_name + "\n"
        " - old size: " + std::to_string(decomp->size) + "\n"
        " - new size: " + std::to_string(my_offsets.size()) + "\n");
  }
}
/* ----------------------------------------------------------------- */
void enddef (const std::string &filename) {
  auto& s = IOSession::instance();

  auto& f = s.get_file(filename);

  f.set_vars_decomps (s.decomps);

  int ierr = PIOc_enddef(f.ncid);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not end define mode on file.\n"
      " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
      " - file name: " + filename + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");

  f.enddef = true;
}

void redef (const std::string &filename) {
  auto& s = IOSession::instance();

  auto& f = s.get_file(filename);

  int ierr = PIOc_redef(f.ncid);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not resume define mode on file.\n"
      " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
      " - file name: " + filename + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");

  f.enddef = false;
}

void free_unused_decomps () {
  auto& s = IOSession::instance();

  std::vector<std::string> to_be_freed;
  for (auto d : s.decomps) {
    if (d.second->customers==0) {
      to_be_freed.push_back(d.first);
    }
  }

  for (auto d : to_be_freed) {
    s.decomps.erase(d);
  }
}

// ---------- File queries
int get_dim_len (const std::string& filename,
                 const std::string& dimname)
{
  const auto& s = IOSession::instance();
  const auto& f = s.get_file(filename);
  const auto& d = f.get_dim(dimname);
  return detail::get_dim_len(f.ncid,d.ncid);
}

// ----------- Read/Write vars/atts -----------

template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const T& attval)
{
  const auto& s = IOSession::instance();
  const auto& f = s.get_file(filename);

  bool do_enddef = f.enddef;
  if (f.enddef) {
    redef (filename);
  }

  int varid, ierr;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = f.get_var(varname).ncid;
  }
  ierr = PIOc_put_att (f.ncid,varid,attname.c_str(),detail::get_nc_type<T>(),
                       detail::get_len(attval),detail::get_data(attval));
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not set attribute.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");

  if (do_enddef) {
    enddef(filename);
  }
}
// ETI
template void set_attribute<int>(const std::string& filename,
                                 const std::string& varname,
                                 const int& attval);
template void set_attribute<float>(const std::string& filename,
                                   const std::string& varname,
                                   const float& attval);
template void set_attribute<double>(const std::string& filename,
                                    const std::string& varname,
                                    const double& attval);
template void set_attribute<std::string>(const std::string& filename,
                                         const std::string& varname,
                                         const std::string& attval);

template<typename T>
void get_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                          T& attval)
{
  const auto& s = IOSession::instance();
  const auto& f = s.get_file(filename);
  int varid,ierr;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = f.get_var(varname).ncid;
  }

  // Check output value has the right type
  EKAT_REQUIRE_MSG (
      type2str<T>()==detail::nctype2str(detail::get_att_type(f.ncid,varid,attname.c_str())),
      "Error! Attempting to retrieve attribute using the wrong type.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - att name : " + attname + "\n"
      " - att type : " + detail::nctype2str(detail::get_att_type(f.ncid,varid,attname.c_str())) + "\n"
      " - type in  : " + type2str<T>() + "\n");

  // Get attribute length, and resize output if needed
  PIO_Offset len;
  ierr = PIOc_inq_attlen(f.ncid, varid, attname.c_str(), &len);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could inquire attribute length.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
  detail::resize(len,attval);

  // Finally extract att value
  ierr = PIOc_get_att (f.ncid,varid,attname.c_str(),detail::get_data_nonconst(attval));
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not get attribute.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}
// ETI
template void get_attribute<int>(const std::string& filename,
                                 const std::string& varname,
                                       int& attval);
template void get_attribute<float>(const std::string& filename,
                                   const std::string& varname,
                                         float& attval);
template void get_attribute<double>(const std::string& filename,
                                    const std::string& varname,
                                          double& attval);
template void get_attribute<std::string>(const std::string& filename,
                                         const std::string& varname,
                                               std::string& attval);

void update_time(const std::string &filename, const double time)
{
  const auto& s = IOSession::instance();
  const auto& f = s.get_file(filename);
  const auto& t = f.get_var("time");

  PIO_Offset start = detail::get_dim_len(f.ncid,t.ncid);
  PIO_Offset count = 1;
  int ierr;

  // Write time variable
  ierr = PIOc_put_vara(f.ncid,t.ncid,&start,&count,&time);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not update time.\n"
      " - file name: " + filename + "\n"
      " - time val : " + std::to_string(time) + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");

  // Update time frame of all vars
  for (auto it : f.vars) {
    auto& v = *it.second;
    ierr = PIOc_advanceframe(f.ncid,v.ncid);
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not advance time frame of variable.\n"
        " - file name: " + filename + "\n"
        " - var name : " + v.name + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");
  }
}

template<typename T>
void read_variable (const std::string &filename,
                    const std::string &varname,
                    const int time_index,
                          T* buf)
{
  const auto& s = IOSession::instance();
  const auto& f = s.get_file(filename);
  const auto& v = f.get_var(varname);

  int ierr, index = 0;
  const bool has_time_dim =f.dims.find("time")!=f.dims.end();
  if (time_index>=0) {
    EKAT_REQUIRE_MSG (has_time_dim,
        "Error! Time index provide, but no time dimension in file.\n"
        " - file name: " + filename + "\n"
        " - var  name: " + varname + "\n");
    index = time_index;
  } else if (has_time_dim) {
    index = get_dim_len(filename,"time")-1;
  }

  if (index>0) {
    ierr = PIOc_setframe(f.ncid, v.ncid, index);

    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could set time frame.\n"
        " - file name: " + filename + "\n"
        " - var name : " + varname + "\n"
        " - t index  : " + std::to_string(time_index) + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");
  }

  ierr = PIOc_read_darray(f.ncid, v.ncid,
                          v.decomp->ncid, v.decomp->size,
                          buf);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could get attribute.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}
// ETI
template void read_variable<int>(const std::string& filename,
                                 const std::string& varname,
                                 const int time_index,
                                       int* buf);
template void read_variable<float>(const std::string& filename,
                                   const std::string& varname,
                                   const int time_index,
                                         float* buf);
template void read_variable<double>(const std::string& filename,
                                    const std::string& varname,
                                    const int time_index,
                                          double* buf);

template<typename T>
void write_variable (const std::string &filename,
                     const std::string &varname,
                     const T* buf)
{
  const auto& s = IOSession::instance();
  const auto& f = s.get_file(filename);
  const auto& v = f.get_var(varname);

  EKAT_REQUIRE_MSG (f.enddef,
      "Error! File is not in data mode. Call enddef first.\n"
      " - file name: " + filename + "\n");

  EKAT_REQUIRE_MSG (v.decomp!=nullptr,
      "Error! No decomposition set for this variable.\n"
      " - file name: " + filename + "\n"
      " - var  name: " + varname + "\n");

  // Note: scorpio uses void* for write buffers, even though
  //       a const void* should conceptually be ok.
  int ierr = PIOc_write_darray(f.ncid, v.ncid,
                               v.decomp->ncid, v.decomp->size,
                               const_cast<T*>(buf), nullptr);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not write variable.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}
// ETI
template void write_variable<int>(const std::string& filename,
                                  const std::string& varname,
                                  const int* buf);
template void write_variable<float>(const std::string& filename,
                                    const std::string& varname,
                                    const float* buf);
template void write_variable<double>(const std::string& filename,
                                     const std::string& varname,
                                     const double* buf);

} // namespace scorpio
} // namespace scream
