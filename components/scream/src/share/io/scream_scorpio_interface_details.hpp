#ifndef SCREAM_SCORPIO_INTERFACE_DETAILS_HPP
#define SCREAM_SCORPIO_INTERFACE_DETAILS_HPP

#include <pio.h>

#include <string>
#include <vector>

namespace scream {

// ========= Convenience wrappers around PIO_inq_XYZ ============ //

namespace detail {

inline std::string get_att_name (int fid, int attnum) {
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

inline std::vector<std::string> get_var_dims (int fid, int varid) {
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

namespace scorpio {
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

} // namespace scorpio
} // namespace scream

#endif // SCREAM_SCORPIO_INTERFACE_DETAILS_HPP
