#ifndef SCREAM_SCORPIO_INTERFACE_DETAILS_HPP
#define SCREAM_SCORPIO_INTERFACE_DETAILS_HPP

#include <pio.h>

#include <string>
#include <vector>

namespace scream {

// ========= Convenience wrappers around PIO_inq_XYZ ============ //

namespace detail {

inline int str2nctype (const std::string& type) {
  if (type=="int") {
    return PIO_INT;
  } else if (type=="float" || type=="single") {
    return PIO_FLOAT;
  } else if (type=="char") {
    return PIO_CHAR;
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

inline std::string nctype2str (const int nctype) {
  std::string name;
  for (std::string n : {"int","char","float","double"}) {
    if (nctype==str2nctype(n)) {
      name = n;
    }
  }
  EKAT_REQUIRE_MSG (name!="",
      "Error! Unrecognized/unsupported nctype: " + std::to_string(nctype) + "\n");
  return name;
}

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
void resize (PIO_Offset /*len*/, T& /*att*/) {}
template<>
inline void resize<std::string> (PIO_Offset len, std::string& att) { att.resize(len); }

template<typename T>
int get_len (const T&) { return 1; }
template<>
inline int get_len<std::string> (const std::string& s) { return s.size(); }

template<typename T>
const void* get_data (const T& att) { return &att; }
template<>
inline const void* get_data<std::string> (const std::string& s) { return s.data(); }

template<typename T>
void* get_data_nonconst (T& att) { return &att; }
template<>
inline void* get_data_nonconst<std::string> (std::string& s) { return &s[0]; }

inline std::string get_att_name (int fid, int varid, int attnum) {
  std::string name;
  name.resize(PIO_MAX_NAME);
  int ierr = PIOc_inq_attname (fid, varid, attnum, &name[0]);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire attribute name.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - att num: " + std::to_string(attnum) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  // Remove all the \0 trailing chars 
  name = name.c_str();
  return name;
}

inline int get_att_type (int fid, int varid, const char* name) {
  int nctype;
  int ierr = PIOc_inq_atttype (fid, varid, name, &nctype);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire attribute type.\n"
      " - file id : " + std::to_string(fid) + "\n"
      " - var id  : " + std::to_string(varid) + "\n"
      " - att name: " + name + "\n"
      " - pio err : " + std::to_string(ierr) + "\n");
  return nctype;
}

inline std::string get_dim_name (int fid, int dimid) {
  std::string name;
  name.resize(PIO_MAX_NAME);
  int ierr;
  ierr = PIOc_inq_dimname (fid, dimid, &name[0]);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire dimension name.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - dim id : " + std::to_string(dimid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  // Remove all the \0 trailing chars 
  name = name.c_str();
  return name;
}

inline int get_dim_len (int fid, int dimid) {
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

inline std::string get_var_name (int fid, int varid) {
  std::string name;
  name.resize(PIO_MAX_NAME);
  int ierr;
  ierr = PIOc_inq_varname (fid, varid, &name[0], PIO_MAX_NAME);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable name.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  // Remove all the \0 trailing chars 
  name = name.c_str();
  return name;
}

inline int get_var_dtype (int fid, int varid) {
  int nctype,ierr;
  ierr = PIOc_inq_vartype(fid, varid, &nctype);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire var data type.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return nctype;
}

inline int get_var_ndims (int fid, int varid) {
  int ndims, ierr;
  ierr = PIOc_inq_varndims (fid, varid, &ndims);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable number of dimensions.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return ndims;
}

inline int get_var_natts (int fid, int varid) {
  int natts, ierr;
  ierr = PIOc_inq_varnatts (fid, varid, &natts);
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable number of attributes.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");
  return natts;
}

inline std::vector<int> get_var_dims (int fid, int varid) {
  const int ndims = get_var_ndims(fid,varid);
  std::vector<int> dimids (ndims);
  int ierr;
  ierr = PIOc_inq_vardimid (fid, varid, dimids.data());
  EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
      "Error! Could not inquire variable dimensions.\n"
      " - file id: " + std::to_string(fid) + "\n"
      " - var id : " + std::to_string(varid) + "\n"
      " - pio err: " + std::to_string(ierr) + "\n");

  return dimids;
}

inline std::vector<std::string> get_var_dim_names (int fid, int varid) {
  const auto dimids = get_var_dims(fid,varid);

  std::vector<std::string> dims;
  for (const auto& id : dimids) {
    dims.push_back(get_dim_name(fid,id));
  }
  return dims;
}

inline std::vector<int> get_var_dim_lens (int fid, int varid) {
  const auto dimids = get_var_dims(fid,varid);

  std::vector<int> dims;
  for (const auto& id : dimids) {
    dims.push_back(get_dim_len(fid,id));
  }
  return dims;
}

} // namespace detail

namespace scorpio {
// =============== INTERNAL DATA TYPES ================== //

// An entity can be a file, a dim, a var, an attribute,...
struct IOEntity {
  std::string name;
  int ncid = -1;
};

struct IODecomp : IOEntity {
  int customers = 0;
  PIO_Offset size = 0;
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
  bool enddef = false;

  // File content
  strptrmap_t<IOEntity>   dims;
  strptrmap_t<IOVar>      vars;

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

  void load_all_entities_from_file () {
    int ierr, ndims, nvars, ngatts;
    ierr = PIOc_inq(ncid, &ndims, &nvars, &ngatts, nullptr);
    EKAT_REQUIRE_MSG (ierr==PIO_NOERR,
        "Error! Could not inquire pio file.\n"
        " - file name: " + name + "\n"
        " - pio err  : " + std::to_string(ierr) + "\n");

    // Add dims
    for (int idim=0; idim<ndims; ++idim) {
      auto dim = std::make_shared<IOEntity>();
      dim->ncid = idim;
      dim->name = detail::get_dim_name(ncid,dim->ncid);
      dims[dim->name] = dim;
    }

    // Add vars
    for (int ivar=0; ivar<nvars; ++ivar) {
      auto var = std::make_shared<IOVar>();
      var->ncid = ivar;
      var->name = detail::get_var_name(ncid,var->ncid);
      vars[var->name] = var;
    }
  }

  void set_vars_decomps (const strptrmap_t<IODecomp>& decomps) {
    // Adjust the number of customers of stored decomps
    for (auto it : vars) {
      auto& v = *it.second;
      auto dim_names = detail::get_var_dim_names (ncid,v.ncid);
      auto dim_lens  = detail::get_var_dim_lens  (ncid,v.ncid);
      if (dim_names[0]=="time") {
        // We always read/write one time slice at a time.
        dim_lens[0] = 1;
      }
      auto dtype = detail::get_var_dtype (ncid,v.ncid);

      auto decomp_tag = detail::nctype2str (dtype) + "."
                      + ekat::join(dim_names,"-") + "."
                      + ekat::join(dim_lens,"-");

      auto d_it = decomps.find(decomp_tag);
      if (d_it!=decomps.end()) {
        auto decomp = d_it->second;
        v.decomp = decomp;
        ++decomp->customers;
      }
    }
  }
};

// Internal singleton to persistently store some PIO stuff
struct IOSession {
  static IOSession& instance (const bool must_be_inited = true) {
    static IOSession s;
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
  IOSession () = default;
};

} // namespace scorpio
} // namespace scream

#endif // SCREAM_SCORPIO_INTERFACE_DETAILS_HPP
