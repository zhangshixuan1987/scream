#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "ekat/util/ekat_string_utils.hpp"
#include "share/field/field_tag.hpp"
#include "share/scream_types.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <vector>

/* C++/F90 bridge to F90 SCORPIO routines */

namespace scream {
namespace scorpio {

using offset_t = std::int64_t;

// Retrieve the int codes PIO uses to specify data types
int str2nctype (const std::string& type);

// WARNING: these values must match the ones of file_purpose_in and file_purpose_out
// in the scream_scorpio_interface F90 module
enum class FileMode {
  Read = 1,
  Write = 2
};

// Initialize/finalize PIO
void init (const ekat::Comm& comm, const int atm_id = 0);
bool is_inited();
void finalize ();

// File handling
void open_file (const std::string& filename, const FileMode mode);
void release_file (const std::string& filename);

void register_dimension (const std::string& filename,
                         const std::string& name,
                         const int length);
void register_variable (const std::string& filename,
                        const std::string& shortname,
                        const std::string& longname,
                        const std::string& units,
                        const std::vector<std::string>& dimensions,
                        const std::string& dtype);

inline
void register_variable (const std::string& filename,
                        const std::string& name,
                        const std::string& units,
                        const std::vector<std::string>& dimensions,
                        const std::string& dtype) {
  register_variable(filename,name,"",units,dimensions,dtype);
}
inline
void register_variable (const std::string& filename,
                        const std::string& name,
                        const std::vector<std::string>& dimensions,
                        const std::string& dtype) {
  register_variable(filename,name,"","",dimensions,dtype);
}

void register_decomp (const std::string& dtype,
                      const std::vector<std::string>& dims,
                      const std::vector<int>& gdimlen,
                      const std::vector<offset_t>& my_offsets);

void enddef (const std::string &filename);
void redef (const std::string &filename);
void free_unused_decomps ();

// File queries
int get_dim_len (const std::string& filename,
                 const std::string& dimname);

// Var/att read/write
template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const T& attval);
template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& attname,
                    const T& attval)
{
  std::string global = "GLOBAL";
  set_attribute(filename,global,attname,attval);
}
// Overload, since string literals are not caught
template<std::size_t N>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const char (&attval)[N]) {
  set_attribute<std::string>(filename,varname,attname,attval);
}
template<std::size_t N>
void set_attribute (const std::string& filename,
                    const std::string& attname,
                    const char (&attval)[N]) {
  set_attribute<std::string>(filename,attname,attval);
}

template<typename T>
void get_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                          T& attval);

template<typename T>
void get_attribute (const std::string& filename,
                    const std::string& attname,
                          T& attval)
{
  std::string global = "GLOBAL";
  get_attribute(filename,global,attname,attval);
}

void update_time(const std::string &filename, const double time);

template<typename T>
void read_variable (const std::string &filename,
                    const std::string &varname,
                    const int time_index,
                          T* buf);
template<typename T>
void read_variable (const std::string &filename,
                    const std::string &varname,
                          T* buf)
{
  read_variable (filename,varname,0,buf);
}

template<typename T>
void write_variable (const std::string &filename,
                     const std::string &varname,
                     const T* buf);

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
