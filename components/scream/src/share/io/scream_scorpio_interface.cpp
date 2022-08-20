#include "scream_scorpio_interface.hpp"

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
namespace detail{
// Query PIO using raw ids
std::string get_att_name (int fid, int attnum);
std::string get_dim_name (int fid, int dimid);
int get_dim_len (int fid, int dimid);
std::string get_var_name (int fid, int varid);
int get_var_dtype (int fid, int varid);
int get_var_ndims (int fid, int varid);
std::vector<std::string> get_var_dims (int fid, int varid);
} // namespace detail

namespace scorpio {

inline std::string e2str (const FileMode mode) {
  return mode==FileMode::Read ? "Read" : "Write";
}

int str2nctype (const std::string& type) {
  if (type=="int") {
    return PIO_INT;
  } else if (type=="float" || type=="single") {
    return PIO_FLOAT;
  } else if (type=="double") {
    return PIO_DOUBLE;
  } else if (type=="real") {
#if defined(SCREAM_DOUBLE_PRECISION)
    return PIO_DOUBLE;
#else
    return PIO_FLOAT;
#endif
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type '" + type + "'.\n");
  }
}

std::string nctype2str (const int nctype) {
  std::string name;
  for (std::string n : {"int","float","double"}) {
    if (nctype==str2nctype(n)) {
      name = n;
    }
  }
  EKAT_REQUIRE_MSG (name!="",
      "Error! Unrecognized/unsupported nctype: " + std::to_string(nctype) + "\n");
  return name;
}

// =============== INTERNAL DATA TYPES ================== //

// An entity can be a file, a dim, a var, an attribute,...
struct IOEntity {
  std::string name;
  int ncid = -1;
};

struct IODecomp : IOEntity {
  int customers = 0;
  int size = 0;
};

struct IOVar : IOEntity {
  std::shared_ptr<IODecomp> decomp;
};

template<typename T>
using strptrmap_t = std::map<std::string,std::shared_ptr<T>>;

// Container for dims/vars/atts
struct IOFile : IOEntity {
  FileMode      mode;
  int customers = 0;

  // File content
  strptrmap_t<IOEntity>   dims;
  strptrmap_t<IOVar>      vars;
  strptrmap_t<IOEntity>   atts;

  const IOEntity& get_dim (const std::string& dname) const {
    EKAT_REQUIRE_MSG (dims.find(dname)!=dims.end(),
        "Error! Could not find dimension in file.\n"
        "   - file name: " + name + "\n"
        "   - dim  name: " + dname + "\n");
    return *dims.at(dname);
  }
  const IOVar& get_var (const std::string& vname) const {
    EKAT_REQUIRE_MSG (vars.find(vname)!=vars.end(),
        "Error! Could not find variable in file.\n"
        "   - file name: " + name + "\n"
        "   - var  name: " + vname + "\n");
    return *vars.at(vname);
  }
  const IOEntity& get_att (const std::string& attname) const {
    EKAT_REQUIRE_MSG (atts.find(attname)!=atts.end(),
        "Error! Could not find attribute in file.\n"
        "   - file name: " + name + "\n"
        "   - att  name: " + attname + "\n");
    return *atts.at(attname);
  }

  void load_all_entities_from_file () {
    int ierr, ndims, nvars, ngatts;
    ierr = PIOc_inq(ncid, &ndims, &nvars, &ngatts, nullptr);
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not inquire pio file.\n"
        " - file name: " + name + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    // Add dims
    for (int i=0; i<ndims; ++i) {
      auto dim = std::make_shared<IOEntity>();
      dim->ncid = i;
      dim->name = detail::get_dim_name(ncid,dim->ncid);
      dims[dim->name] = dim;
    }

    // Add vars
    for (int i=0; i<nvars; ++i) {
      auto var = std::make_shared<IOVar>();
      var->ncid = i;
      var->name = detail::get_var_name(ncid,var->ncid);
      vars[var->name] = var;
    }

    // Add atts
    for (int i=0; i<ngatts; ++i) {
      auto att = std::make_shared<IOEntity>();
      att->ncid = PIO_GLOBAL;
      att->name = detail::get_att_name(ncid,i);
      atts[att->name] = att;
    }
  }
};

// Internal singleton to persistently store some PIO stuff
struct ScorpioSession {
  static ScorpioSession& instance (const bool must_be_inited = true) {
    static ScorpioSession s;
    if (must_be_inited) {
      EKAT_REQUIRE_MSG (s.iosysid>=0,
          "Error! PIO subsystem was expected to be already inited.\n");
    } else {
      EKAT_REQUIRE_MSG (s.iosysid==-1,
          "Error! PIO subsystem was expected to not be inited yet.\n");
    }
    return s;
  }

  int iosysid = -1;
  int iotype  = -1;
  int rearr   = -1;
  int format  = -1;

  bool must_finalize = false;

  const IOFile& get_file (const std::string& name) const {
    EKAT_REQUIRE_MSG (files.find(name)!=files.end(),
        "Error! Could not find file '" + name + "'.\n");
    return *files.at(name);
  }

  IOFile& get_file (const std::string& name) {
    EKAT_REQUIRE_MSG (files.find(name)!=files.end(),
        "Error! Could not find file '" + name + "'.\n");
    return *files.at(name);
  }

  strptrmap_t<IOFile>     files;
  strptrmap_t<IODecomp>   decomps;

  ekat::Comm comm;

  void remove_file (const std::string fname);

private:
  ScorpioSession () = default;
};

// ============================================================== //

// ------------------------ IMPLEMENTATION ------------------------- //

void init(const ekat::Comm& comm, const int atm_id) {
  auto& s = ScorpioSession::instance(false);

  s.comm = comm;
#ifdef SCREAM_CIME_BUILD
  // PIO system should alreay be inited. Grab it.
  s.iosysid = get_io_sys(atm_id);
  s.iotype  = get_io_type(atm_id);
  s.rearr   = get_io_rearranger(atm_id);
  s.format  = get_io_format(atm_id);
#else
  s.iotype = PIO_IOTYPE_PNETCDF;
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
  auto& s = ScorpioSession::instance();

  if (s.must_finalize) {
    int err = PIOc_finalize(s.iosysid);

    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Could not finalize pio subsystem.\n"
        " - atm rank: " + std::to_string(s.comm.rank()) + "\n"
        " - pio err : " + std::to_string(err) + "\n");
  }
}

void register_file (const std::string& filename, const FileMode mode) {
  auto& s = ScorpioSession::instance();

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
    } else {
      ierr = PIOc_createfile(s.iosysid,&f->ncid,&s.iotype,filename.c_str(),PIO_WRITE | PIO_CLOBBER);
    }
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not open pio file.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");
  }

  ++f->customers;
}
/* ----------------------------------------------------------------- */
void closefile (const std::string& filename) {
  auto& s = ScorpioSession::instance();

  EKAT_REQUIRE_MSG (s.files.find(filename)!=s.files.end(),
      "Error! Could not close file, since it was not found.\n"
      " - file name: " + filename + "\n");

  auto& f = s.get_file(filename);
  --f.customers;

  if (f.customers==0) {
    // Remove vars as customers of their decomp
    for (auto v : f.vars) {
      --v.second->decomp->customers;
      // auto dims  = detail::get_var_dims  (f.ncid,v.second->ncid);
      // auto dtype = detail::get_var_dtype (f.ncid,v.second->ncid);

      // auto decomp_tag = nctype2str (dtype) + "." + ekat::join(dims,"-");
      // EKAT_REQUIRE_MSG (s.decomps.find(decomp_tag)!=s.decomps.end(),
      //     "Error! Decomp tag for a registered variable was already removed.\n"
      //     " - atm rank  : " + std::to_string(s.comm.rank()) + "\n"
      //     " - file name : " + filename + "\n"
      //     " - var name  : " + v.second->name + "\n"
      //     " - decomp tag: " + decomp_tag + "\n");

      // auto decomp = s.decomps.at(decomp_tag);
      // --decomp->customers;
    }

    // File is no longer needed. Clean it up
    s.files.erase(filename);
  }
}
/* ------------------------------------------------- */
void register_dimension (const std::string& filename,
                         const std::string& dim_name,
                         const int length)
{
  auto& s = ScorpioSession::instance();

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
                        const std::string& var_name,
                        const std::string& longname,
                        const std::string& units,
                        const std::vector<std::string>& dim_names,
                        const std::string& dtype)
{
  auto& s = ScorpioSession::instance();

  int ierr;
  auto& f = s.get_file(filename);
  auto& v = f.vars[var_name];
  if (v->ncid == -1) {
    EKAT_REQUIRE_MSG (f.mode==FileMode::Write,
        "Error! Registering new variable in a read-only file!\n");

    std::vector<int> dimids;
    for (const auto& dn : dim_names) {
      dimids.push_back(f.get_dim(dn).ncid);
    }
    ierr = PIOc_def_var(f.ncid, var_name.c_str(), str2nctype(dtype),
                        dim_names.size(), dimids.data(), &v->ncid);

    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not define variable.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - var name : " + var_name + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    // Add units/longname as attribute
    using vosp = std::vector<std::pair<std::string,std::string>>;
    for (auto it : vosp{ {"units",units}, {"long_name",longname} }) {
      const auto att_n = it.first.c_str();
      const auto att_v = it.second;
      ierr = PIOc_put_att_text(f.ncid, v->ncid, att_n, att_v.size(), att_v.c_str());
      EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
          "Error! Could not add variable attribute.\n"
          " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
          " - file name: " + filename + "\n"
          " - var name : " + var_name + "\n"
          " - attribute: " + att_n + "\n"
          " - pio err  : " + std::to_string(ierr) + "\n");
    }
  } else {
    // Check that var specs are the same
    EKAT_REQUIRE_MSG (detail::get_var_name(f.ncid,v->ncid)==var_name,
        "Error! Something is amiss with dimensions storage.\n"
        "       Please, contact developers.\n");
    const auto nctype = detail::get_var_dtype(f.ncid,v->ncid);
    EKAT_REQUIRE_MSG (nctype==str2nctype(dtype),
        "Error! Variable already defined with different data type.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - var name : " + var_name + "\n"
        " - var type : " + nctype2str(nctype) + "\n"
        " - type in  : " + dtype + "\n");
    EKAT_REQUIRE_MSG (dim_names==detail::get_var_dims(f.ncid,v->ncid),
        "Error! Variable already defined with different dimensions.\n"
        " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
        " - file name: " + filename + "\n"
        " - var name : " + var_name + "\n"
        " - var dims : (" + ekat::join(detail::get_var_dims(f.ncid,v->ncid),", ") + ")\n"
        " - dims in  : (" + ekat::join(dim_names,", ") + ")\n");
  }
}
/* ----------------------------------------------------------------- */
void register_decomp (const std::string& dtype,
                      const std::vector<std::string>& dims,
                      const std::vector<int>& gdimlen,
                      const std::vector<offset_t>& my_offsets)
{
  auto& s = ScorpioSession::instance();

  EKAT_REQUIRE_MSG (dims.size()==gdimlen.size(),
      "[register_decomp] Error! Input arrays have different lengths.\n"
      " - dims size   : " + std::to_string(dims.size()) + "\n"
      " - gdimlen size: " + std::to_string(gdimlen.size()) + "\n");

  int ierr;
  const auto decomp_name = dtype + "." + ekat::join(dims,"-");

  auto& decomp = s.decomps[decomp_name];
  if (decomp==nullptr) {
    // Create decomp
    decomp = std::make_shared<IODecomp>();
    decomp->name = decomp_name;
    decomp->size = my_offsets.size();

    auto offsets = reinterpret_cast<const PIO_Offset*>(my_offsets.data());
    ierr = PIOc_InitDecomp(s.iosysid, str2nctype(dtype), dims.size(),
                           gdimlen.data(), my_offsets.size(),
                           offsets, &decomp->ncid, &s.rearr,
                           nullptr, nullptr);

    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not add decomposition.\n"
        " - atm rank: " + std::to_string(s.comm.rank()) + "\n"
        " - decomp  : " + decomp_name + "\n"
        " - pio err : " + std::to_string(ierr) + "\n");
  }
}
/* ----------------------------------------------------------------- */
void enddef (const std::string &filename) {
  auto& s = ScorpioSession::instance();

  auto& f = s.get_file(filename);

  // Adjust the number of customers of stored decomps
  for (auto it : f.vars) {
    auto& v = *it.second;
    auto dims  = detail::get_var_dims  (f.ncid,v.ncid);
    auto dtype = detail::get_var_dtype (f.ncid,v.ncid);

    auto decomp_tag = nctype2str (dtype) + "." + ekat::join(dims,"-");
    EKAT_REQUIRE_MSG (s.decomps.find(decomp_tag)!=s.decomps.end(),
        "Error! Missing decomp tag for a registered variable.\n"
        " - atm rank  : " + std::to_string(s.comm.rank()) + "\n"
        " - file name : " + filename + "\n"
        " - var name  : " + v.name + "\n"
        " - decomp tag: " + decomp_tag + "\n");

    v.decomp = s.decomps.at(decomp_tag); 
    ++v.decomp->customers;
  }

  int ierr = PIOc_enddef(f.ncid);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not end define mode on file.\n"
      " - atm rank : " + std::to_string(s.comm.rank()) + "\n"
      " - file name: " + filename + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}

void free_unused_decomps () {
  auto& s = ScorpioSession::instance();

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
  const auto& s = ScorpioSession::instance();
  const auto& f = s.get_file(filename);
  const auto& d = f.get_dim(dimname);
  return detail::get_dim_len(f.ncid,d.ncid);
}

bool is_file_open (const std::string& filename) {
  const auto& s = ScorpioSession::instance();
  return s.files.find(filename)!=s.files.end();
}

// ----------- Read/Write vars/atts -----------

template<typename T>
int get_nc_type () {
  int nctype;
  if (std::is_same<double,T>::value) {
    nctype = PIO_DOUBLE;
  } else if (std::is_same<int,T>::value) {
    nctype = PIO_INT;
  } else if (std::is_same<float,T>::value) {
    nctype = PIO_FLOAT;
  } else if (std::is_same<std::string,T>::value) {
    nctype = PIO_CHAR;
  } else {
    EKAT_ERROR_MSG ("Error! Unsupported/unknown type.\n");
  }
  return nctype;
}

template<typename T>
void resize_att (PIO_Offset /*len*/, T& /*att*/) {
  // Nothing to do
}

template<>
void resize_att<std::string> (PIO_Offset len, std::string& att) {
  att.resize(len);
}

template<typename T>
int get_att_len (const T&) {
  return 1;
}
template<>
int get_att_len<std::string> (const std::string& s) {
  return s.size();
}
template<typename T>
const void* get_att_data (const T& att) {
  return &att;
}
template<>
const void* get_att_data<std::string> (const std::string& s) {
  return s.data();
}
template<typename T>
void* get_att_data_nonconst (T& att) {
  return &att;
}
template<>
void* get_att_data_nonconst<std::string> (std::string& s) {
  return &s[0];
}

template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& att_name,
                    const T& att_val)
{
  const auto& s = ScorpioSession::instance();
  const auto& f = s.get_file(filename);
  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = f.get_var(varname).ncid;
  }
  int ierr = PIOc_put_att (f.ncid,varid,att_name.c_str(),get_nc_type<T>(),
                           get_att_len(att_val),get_att_data(att_val));
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not set attribute.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}
// ETI
template void set_attribute<int>(const std::string& filename,
                                 const std::string& varname,
                                 const int& att_val);
template void set_attribute<float>(const std::string& filename,
                                   const std::string& varname,
                                   const float& att_val);
template void set_attribute<double>(const std::string& filename,
                                    const std::string& varname,
                                    const double& att_val);
template void set_attribute<std::string>(const std::string& filename,
                                         const std::string& varname,
                                         const std::string& att_val);

template<typename T>
void get_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& att_name,
                          T& att_val)
{
  const auto& s = ScorpioSession::instance();
  const auto& f = s.get_file(filename);
  int varid,ierr;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = f.get_var(varname).ncid;
  }
  PIO_Offset len;
  ierr = PIOc_inq_attlen(f.ncid, varid, att_name.c_str(), &len);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could inquire attribute length.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");

  resize_att(len,att_val);

  ierr = PIOc_get_att (f.ncid,varid,att_name.c_str(),get_att_data_nonconst(att_val));
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not get attribute.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}
// ETI
template void get_attribute<int>(const std::string& filename,
                                 const std::string& varname,
                                       int& att_val);
template void get_attribute<float>(const std::string& filename,
                                   const std::string& varname,
                                         float& att_val);
template void get_attribute<double>(const std::string& filename,
                                    const std::string& varname,
                                          double& att_val);
template void get_attribute<std::string>(const std::string& filename,
                                         const std::string& varname,
                                               std::string& att_val);

void update_time(const std::string &filename, const double time)
{
  const auto& s = ScorpioSession::instance();
  const auto& f = s.get_file(filename);
  const auto& t = f.get_var("time");

  int ierr = PIOc_put_var(f.ncid,t.ncid,&time);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not update time.\n"
      " - file name: " + filename + "\n"
      " - time val : " + std::to_string(time) + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");
}

template<typename T>
void read_variable (const std::string &filename,
                    const std::string &varname,
                    const int time_index,
                          T* buf)
{
  const auto& s = ScorpioSession::instance();
  const auto& f = s.get_file(filename);
  const auto& v = f.get_var(varname);

  int ierr, index;
  if (time_index>0) {
    index = time_index;
  } else {
    index = get_dim_len(filename,"time");
  }

  ierr = PIOc_setframe(f.ncid, v.ncid, index);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could set time frame.\n"
      " - file name: " + filename + "\n"
      " - var name : " + varname + "\n"
      " - t index  : " + std::to_string(time_index) + "\n"
      " - pio err  : " + std::to_string(ierr) + "\n");

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
  const auto& s = ScorpioSession::instance();
  const auto& f = s.get_file(filename);
  const auto& v = f.get_var(varname);

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

// ========= Convenience wrappers around PIO_inq_XYZ ============ //

namespace detail {
std::string get_att_name (int fid, int attnum) {
  std::string name;
  name.resize(PIO_MAX_NAME);
  int ierr = PIOc_inq_attname (fid, PIO_GLOBAL, attnum, &name[0]);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire global attribute name.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - att num: " + std::to_string(attnum) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return name;
}


std::string get_dim_name (int fid, int dimid) {
  std::string name;
  name.resize(PIO_MAX_NAME);
  int ierr;
  ierr = PIOc_inq_dimname (fid, dimid, &name[0]);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire dimension name.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - dim id : " + std::to_string(dimid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return name;
}

int get_dim_len (int fid, int dimid) {
  PIO_Offset len;
  int ierr;
  ierr = PIOc_inq_dimlen (fid, dimid, &len);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire dimension length.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - dim id : " + std::to_string(dimid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return len;
}

std::string get_var_name (int fid, int varid) {
  std::string name;
  name.resize(PIO_MAX_NAME);
  int ierr;
  ierr = PIOc_inq_varname (fid, varid, &name[0], PIO_MAX_NAME);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable name.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return name;
}

int get_var_dtype (int fid, int varid) {
  int nctype,ierr;
  ierr = PIOc_inq_vartype(fid, varid, &nctype);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire var data type.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return nctype;
}

int get_var_ndims (int fid, int varid) {
  int ndims, ierr;
  ierr = PIOc_inq_varndims (fid, varid, &ndims);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable number of dimensions.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return ndims;
}

std::vector<std::string> get_var_dims (int fid, int varid) {
  const int ndims = detail::get_var_ndims(fid,varid);
  std::vector<int> dimids (ndims);
  int ierr;
  ierr = PIOc_inq_vardimid (fid, varid, dimids.data());
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable dimensions.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");

  std::vector<std::string> dims;
  for (const auto& id : dimids) {
    dims.push_back(detail::get_dim_name(fid,id));
  }
  return dims;
}

} // namespace detail
} // namespace scream
