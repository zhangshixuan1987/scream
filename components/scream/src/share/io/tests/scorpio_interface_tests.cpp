#include <catch2/catch.hpp>

#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scream_scorpio_interface_details.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include <ekat/util/ekat_test_utils.hpp>

using namespace scream;

TEST_CASE("file-manip")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  namespace S = scorpio;
  using FM = S::FileMode;

  auto engine = setup_random_test ();
  using IPDF = std::uniform_int_distribution<int>;
  using FPDF = std::uniform_real_distribution<float>;
  using DPDF = std::uniform_real_distribution<double>;

  S::init(comm);

  // Make sure we can test all internal errors,
  // without scorpio just killing the execution.
  auto&s = S::IOSession::instance();
  PIOc_Set_IOSystem_Error_Handling(s.iosysid, PIO_RETURN_ERROR);

  SECTION ("file") {
    std::string fname = "open_close.nc";

    S::open_file(fname,FM::Write);
    // Cannot have two writes
    REQUIRE_THROWS(S::open_file(fname,FM::Write));
    // Cannot read/write at the same time
    REQUIRE_THROWS(S::open_file(fname,FM::Read));
    S::release_file(fname);

    // Double read is ok. Open/close twice, to emulate
    // two objects separatedly reading
    S::open_file(fname,FM::Read);
    S::open_file(fname,FM::Read);
    S::release_file(fname);
    S::release_file(fname);
    // Last release_file has closed the file
    REQUIRE_THROWS(S::release_file(fname));
  }

  SECTION ("register") {
    std::string fname = "register.nc";
    S::open_file(fname,FM::Write);

    // Dimensions
    S::register_dimension(fname,"dim1",4);
    S::register_dimension(fname,"dim2",3);
    S::register_dimension(fname,"dim1",4);
    // Only one unlimited dimension allowed
    REQUIRE_THROWS(S::register_dimension(fname,"time2",PIO_UNLIMITED));
    // Already registered with different dim
    REQUIRE_THROWS (S::register_dimension(fname,"dim1",12));

    // Variables
    S::register_variable(fname,"var1",{"dim1","dim2"},"int");
    S::register_variable(fname,"var2",{"dim1","dim2"},"double");
    S::register_variable(fname,"var2",{"dim1","dim2"},"double");
    S::register_variable(fname,"var3",{"dim1","dim2"},"float");
    // Already registered, with different data type
    REQUIRE_THROWS (S::register_variable(fname,"var1",{"dim1","dim2"},"double"));
    // Already registered, with different dims
    REQUIRE_THROWS (S::register_variable(fname,"var1",{"dim1","dim1"},"int"));
    // Unsupported data type
    REQUIRE_THROWS (S::register_variable(fname,"var3",{"dim1","dim1"},"char"));
    // Unrecognized dim
    REQUIRE_THROWS (S::register_variable(fname,"var3",{"dim1","dim3"},"char"));

    S::release_file(fname);
  }

  SECTION ("variables") {
    std::string fname = "vars.nc";
    S::open_file(fname,FM::Write);

    const int dim1 = 4;
    const int dim2 = 3;
    const int size = dim1*dim2;

    const int gdim1 = dim1*comm.size();
    std::vector<S::offset_t> my_offsets(size);
    std::iota(my_offsets.begin(),my_offsets.end(),size*comm.rank());

    // Register
    S::register_dimension(fname,"dim1",gdim1);
    S::register_dimension(fname,"dim2",dim2);
    S::register_variable(fname,"ivar",{"dim1","dim2"},"int");
    S::register_variable(fname,"dvar",{"dim1","dim2"},"double");
    S::register_variable(fname,"fvar",{"dim1","dim2"},"float");

    // Generate data
    std::vector<int>    ivar (size);
    std::vector<float>  fvar (size);
    std::vector<double> dvar (size);

    ekat::genRandArray (ivar.data(),size,engine,IPDF(0,100));
    ekat::genRandArray (fvar.data(),size,engine,FPDF(0,100));
    ekat::genRandArray (dvar.data(),size,engine,DPDF(0,100));

    // Set decompositions
    S::register_decomp ("int",{"time","dim1","dim2"},{1,gdim1,dim2},my_offsets);
    S::register_decomp ("float",{"time","dim1","dim2"},{1,gdim1,dim2},my_offsets);
    S::register_decomp ("double",{"time","dim1","dim2"},{1,gdim1,dim2},my_offsets);

    // Not in data mode yet
    REQUIRE_THROWS(S::write_variable(fname,"ivar",ivar.data()));
    S::enddef (fname);

    // Write
    S::write_variable(fname,"ivar",ivar.data());
    S::write_variable(fname,"fvar",fvar.data());
    S::write_variable(fname,"dvar",dvar.data());
    S::release_file(fname);

    // Re-open file
    S::open_file(fname,FM::Read);

    // Read
    std::vector<int>    ivar2 (size,-1);
    std::vector<float>  fvar2 (size,-1);
    std::vector<double> dvar2 (size,-1);
  
    S::read_variable(fname,"ivar",ivar2.data());
    S::read_variable(fname,"fvar",fvar2.data());
    S::read_variable(fname,"dvar",dvar2.data());

    // Check
    for (int i=0; i<size; ++i) {
      REQUIRE (ivar[i]==ivar2[i]);
      REQUIRE (fvar[i]==fvar2[i]);
      REQUIRE (dvar[i]==dvar2[i]);
    }
    S::release_file (fname);
  }

  SECTION ("attributes") {
    std::string fname = "atts.nc";
    S::open_file(fname,FM::Write);

    const int dim1 = 4;
    const int dim2 = 3;
    const int gdim1 = dim1*comm.size();

    // Register
    S::register_dimension(fname,"dim1",gdim1);
    S::register_dimension(fname,"dim2",dim2);
    S::register_variable(fname,"ivar",{"dim1","dim2"},"int");

    // Write atts
    S::set_attribute(fname,"ivar","units","meters");
    S::set_attribute(fname,"author","Luca");
    S::set_attribute(fname,"my_int",42);
    S::enddef (fname);
    S::release_file(fname);

    // Re-open file
    S::open_file(fname,FM::Read);

    // Read atts
    std::string author, units;
    int my_int;
    S::get_attribute(fname,"author",author);
    S::get_attribute(fname,"my_int",my_int);
    S::get_attribute(fname,"ivar","units",units);
    REQUIRE (author=="Luca");
    REQUIRE (my_int==42);
    REQUIRE (units=="meters");

    // Wrong type for the attribute
    REQUIRE_THROWS (S::get_attribute(fname,"my_int",author));
  }
  S::finalize();
}

