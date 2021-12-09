#include <catch2/catch.hpp>
#include <numeric>

#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_property_checks/field_positivity_check.hpp"
#include "share/field/field_property_checks/field_within_interval_check.hpp"
#include "share/field/field_property_checks/field_lower_bound_check.hpp"
#include "share/field/field_property_checks/field_upper_bound_check.hpp"
#include "share/field/field_property_checks/field_nan_check.hpp"
#include "share/field/field_property_checks/check_and_repair_wrapper.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("field_layout") {
  using namespace scream;
  using namespace ShortFieldTagsNames;

  FieldLayout l({EL,GP,GP});

  // Should not be able to set a dimensions vector of wrong rank
  REQUIRE_THROWS(l.set_dimensions({1,2}));

  l.set_dimensions({1,2,3});

  // Should not be able to reset the dimensions once they are set
  REQUIRE_THROWS(l.set_dimensions({1,2,3}));
}

TEST_CASE("field_identifier", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::vector<FieldTag> tags1 = {EL, LEV, CMP};
  std::vector<FieldTag> tags2 = {EL, CMP, LEV};

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 4, 3};

  FieldIdentifier fid1 ("field_1", {tags1,dims1}, kg, "some_grid");
  FieldIdentifier fid2 ("field_1", {tags1,dims1}, kg, "some_grid");
  FieldIdentifier fid3 ("field_1", {tags1,dims2}, kg, "some_grid");
  FieldIdentifier fid4 ("field_2", {tags1,dims2}, kg, "some_grid");
  FieldIdentifier fid5 ("field_2", {tags2,dims2}, kg, "some_grid");
  FieldIdentifier fid6 ("field_2", {tags2,dims2}, m, "some_grid");
  FieldIdentifier fid7 ("field_2", {tags2,dims2}, m, "some_other_grid");

  REQUIRE (fid1==fid2);
  REQUIRE (fid2!=fid3);
  REQUIRE (fid3!=fid4);
  REQUIRE (fid4!=fid5);
  REQUIRE (fid5!=fid6);
  REQUIRE (fid6!=fid7);

  // Check that has_tag option works
  REQUIRE(fid1.get_layout().has_tag(CMP));
  REQUIRE(!fid1.get_layout().has_tag(GP));
}

TEST_CASE("field_tracking", "") {
  using namespace scream;

  FieldTracking track;
  util::TimeStamp time1(2021,10,12,17,8,10);
  util::TimeStamp time2(2021,10,12,17,8,20);
  REQUIRE_NOTHROW (track.update_time_stamp(time2));

  // Cannot rewind time (yet)
  REQUIRE_THROWS  (track.update_time_stamp(time1));
}

TEST_CASE("field", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using kt = KokkosTypes<DefaultDevice>;

  using P4 = ekat::Pack<Real,4>;
  using P8 = ekat::Pack<Real,8>;
  using P16 = ekat::Pack<Real,16>;

  auto engine = setup_random_test ();
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.01,0.99);

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");

  // Check if we can extract a reshaped view
  SECTION ("reshape") {
    Field f1 (fid);

    // Should not be able to reshape before allocating
    REQUIRE_THROWS(f1.get_view<Real*>());

    f1.allocate_view();

    // Reshape should work with both dynamic and static dims
    auto v1 = f1.get_view<Real[3][24]>();
    auto v2 = f1.get_view<Real**>();

    REQUIRE(v1.size()==v2.size());

    // But if wrong static length is used, we should throw
    REQUIRE_THROWS(f1.get_view<Real[3][16]>());

    // Should not be able to reshape to this data type...
    REQUIRE_THROWS(f1.get_view<P16**>());
    // But this should work
    REQUIRE_NOTHROW(f1.get_view<P8**>());

    // Using packs (of allowable size) of different pack sizes
    // should lead to views with different extents.
    // Since there's no padding, their extent on last dimension
    // should be the phys dim divided by pack size.
    auto v3 = f1.get_view<P8**>();
    auto v4 = f1.get_view<P4**>();
    REQUIRE (v4.size() == 2*v3.size());
    REQUIRE (v4.extent_int(0) == fid.get_layout().dim(0));
    REQUIRE (v3.extent_int(0) == fid.get_layout().dim(0));
    REQUIRE (v4.extent_int(1) == fid.get_layout().dim(1) / P4::n);
    REQUIRE (v3.extent_int(1) == fid.get_layout().dim(1) / P8::n);

    // The memory spans should be identical
    REQUIRE (v3.impl_map().memory_span()==v4.impl_map().memory_span());

    // Trying to reshape into something that the allocation cannot accommodate should throw
    REQUIRE_THROWS (f1.get_view<P16***>());
  }

  SECTION ("compare") {

    Field f1(fid), f2(fid);
    f2.get_header().get_alloc_properties().request_allocation(P16::n);
    f1.allocate_view();
    f2.allocate_view();

    auto v1 = f1.get_view<Real**>();
    auto v2 = f2.get_view<P8**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = i*dim1+j;

      int jpack = j / P8::n;
      int jvec = j % P8::n;
      v2(i,jpack)[jvec] = i*dim1+j;
    });
    Kokkos::fence();

    // The views were filled the same way, so they should test equal
    // NOTE: this cmp function only test the "actual" field, discarding padding.
    REQUIRE(views_are_equal(f1,f2));

    // Check self equivalence
    // get_const returns a copy of self, so equivalent (if already allocated)
    REQUIRE (f1.equivalent(f1.get_const()));
    REQUIRE (f1.equivalent(f1));
    // f1 and f2 have independent views, so they are not equivalent.
    REQUIRE (!f1.equivalent(f2));
  }

  // Check copy constructor
  SECTION ("copy ctor") {
    Field f1 (fid);

    f1.allocate_view();
    f1.deep_copy(3.0);

    Field f2 = f1;
    REQUIRE(f2.get_header_ptr()==f1.get_header_ptr());
    REQUIRE(f2.get_internal_view_data<Real>()==f1.get_internal_view_data<Real>());
    REQUIRE(f2.is_allocated());
    REQUIRE(views_are_equal(f1,f2));
  }

  SECTION ("deep_copy") {
    std::vector<FieldTag> t1 = {COL,CMP,LEV};
    std::vector<int> d1 = {3,2,24};

    FieldIdentifier fid1("vec_3d",{t1,d1},m/s,"some_grid");

    Field f1(fid1);
    f1.allocate_view();
    f1.deep_copy(1.0);
    f1.sync_to_host();
    auto v = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    for (int i=0; i<fid1.get_layout().size(); ++i) {
      REQUIRE (v[i]==1.0);
    }
  }

  // Subfields
  SECTION ("subfield") {
    std::vector<FieldTag> t1 = {COL,CMP,CMP,LEV};
    std::vector<int> d1 = {3,10,2,24};

    FieldIdentifier fid1("4d",{t1,d1},m/s,"some_grid");

    Field f1(fid1);
    f1.allocate_view();
    randomize(f1,engine,pdf);

    const int idim = 1;
    const int ivar = 2;

    auto f2 = f1.subfield(idim,ivar);

    // Wrong rank for the subfield f2
    REQUIRE_THROWS(f2.get_view<Real****>());

    auto v4d_h = f1.get_view<Real****,Host>();
    auto v3d_h = f2.get_view<Real***,Host>();
    for (int i=0; i<d1[0]; ++i)
      for (int j=0; j<d1[2]; ++j)
        for (int k=0; k<d1[3]; ++k) {
          REQUIRE (v4d_h(i,ivar,j,k)==v3d_h(i,j,k));
        }
  }

  // Dynamic Subfields
  SECTION ("dynamic_subfield") {
    const int vec_dim = 10;
    std::vector<FieldTag> t1 = {COL,CMP,CMP,LEV};
    std::vector<int> d1 = {3,vec_dim,2,24};

    FieldIdentifier fid1("4d",{t1,d1},m/s,"some_grid");

    Field f1(fid1);
    f1.allocate_view();
    randomize(f1,engine,pdf);

    const int idim = 1;
    const int ivar = 0;

    auto f2 = f1.subfield(idim,ivar,/* dynamic = */ true);

    // Cannot reset subview idx of non-subfield fields
    REQUIRE_THROWS(f1.get_header().get_alloc_properties().reset_subview_idx(0));

    // subview idx out of bounds
    auto& f2_ap = f2.get_header().get_alloc_properties();
    REQUIRE_THROWS(f2_ap.reset_subview_idx(-1));
    REQUIRE_THROWS(f2_ap.reset_subview_idx(vec_dim));

    // Fill f1 with random numbers, and verify corresponding subviews get same values
    randomize(f1,engine,pdf);

    for (int ivar_dyn=0; ivar_dyn<vec_dim; ++ivar_dyn) {
      // Reset slice idx
      f2_ap.reset_subview_idx(ivar_dyn);
      REQUIRE(f2_ap.get_subview_info().slice_idx==ivar_dyn);

      auto v4d_h = f1.get_view<Real****,Host>();
      auto v3d_h = f2.get_view<Real***,Host>();
      for (int i=0; i<d1[0]; ++i)
        for (int j=0; j<d1[2]; ++j)
          for (int k=0; k<d1[3]; ++k) {
            REQUIRE (v4d_h(i,ivar_dyn,j,k)==v3d_h(i,j,k));
          }
    }
  }

  SECTION ("vector_component") {
    std::vector<FieldTag> tags_2 = {COL,CMP,LEV};
    std::vector<int> dims_2 = {3,2,24};

    FieldIdentifier fid_2("vec_3d",{tags_2,dims_2},m/s,"some_grid");

    Field f_vec(fid_2);
    f_vec.allocate_view();

    auto f0 = f_vec.get_component(0);
    auto f1 = f_vec.get_component(1);

    // No 3rd component
    REQUIRE_THROWS(f_vec.get_component(2));

    // f0 is scalar, no vector dimension
    REQUIRE_THROWS(f0.get_component(0));

    f0.deep_copy(1.0);
    f1.deep_copy(2.0);

    f_vec.sync_to_host();

    auto v = f_vec.get_view<Real***,Host>();
    for (int col=0; col<3; ++col) {
      for (int lev=0; lev<24; ++lev) {
        REQUIRE (v(col,0,lev)==1.0);
        REQUIRE (v(col,1,lev)==2.0);
      }
    }
  }


  SECTION ("host_view") {
    Field f(fid);

    // Views not yet allocated
    REQUIRE_THROWS(f.get_internal_view_data<Real>());
    REQUIRE_THROWS(f.get_internal_view_data<Real,Host>());
    REQUIRE_THROWS(f.sync_to_host());
    REQUIRE_THROWS(f.sync_to_dev());

    f.allocate_view();
    randomize(f,engine,pdf);

    // Get reshaped view on device, and manually create Host mirror
    auto v2d = f.get_view<Real**>();
    auto v2d_hm = Kokkos::create_mirror_view(v2d);
    Kokkos::deep_copy(v2d_hm,v2d);

    // Get reshaped view straight on Host
    auto v2dh = f.get_view<Real**,Host>();

    // The two should match
    for (int i=0; i<dims[0]; ++i) {
      for (int j=0; j<dims[1]; ++j) {
        REQUIRE (v2dh(i,j) == v2d_hm(i,j) );
      }
    }
  }
}

TEST_CASE("field_property_check", "") {

  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::vector<FieldTag> tags = {EL, GP, LEV};
  std::vector<int> dims = {2, 3, 12};

  FieldIdentifier fid ("field_1",{tags,dims}, m/s,"some_grid");

  auto engine = setup_random_test();
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pos_pdf(0.01,0.99);
  RPDF neg_pdf(-0.99, -0.01);

  // Check positivity.
  SECTION ("field_positivity_check") {
    Field f1(fid);
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    // Assign positive values to the field and make sure it passes our test for
    // positivity.
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f1_data,num_reals,engine,pos_pdf);
    f1.sync_to_dev();

    auto positivity_check = std::make_shared<FieldPositivityCheck>(false);
    REQUIRE(not positivity_check->can_repair());
    REQUIRE(positivity_check->check(f1));

    // Assign non-positive values to the field and make sure it fails the check.
    ekat::genRandArray(f1_data,num_reals,engine,neg_pdf);
    f1.sync_to_dev();
    REQUIRE(not positivity_check->check(f1));
  }

  // Check positivity with repairs.
  SECTION ("field_positivity_check_with_repairs") {
    Field f1(fid);
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    // Assign non-positive values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f1_data,num_reals,engine,neg_pdf);
    f1.sync_to_dev();

    auto positivity_check = std::make_shared<FieldPositivityCheck>(1);
    REQUIRE(positivity_check->can_repair());
    REQUIRE(not positivity_check->check(f1));
    positivity_check->repair(f1);
    REQUIRE(positivity_check->check(f1));
  }

  // Check that values are not NaN
  SECTION("field_not_nan_check") {
    Field f1(fid);
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    auto nan_check = std::make_shared<FieldNaNCheck>();

    // Assign  values to the field and make sure it passes our test for NaNs.
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f1_data,num_reals,engine,neg_pdf);
    f1.sync_to_dev();
    REQUIRE(nan_check->check(f1));

    // Assign a NaN value to the field, make sure it fails the check,
    Int midpt = num_reals / 2;
    f1_data[midpt] = std::numeric_limits<Real>::quiet_NaN();
    f1.sync_to_dev();
    REQUIRE(not nan_check->check(f1));
  }

  // Check that the values of a field lie within an interval.
  SECTION ("field_within_interval_check") {
    Field f1(fid);
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    auto interval_check = std::make_shared<FieldWithinIntervalCheck>(0, 1);
    REQUIRE(interval_check->can_repair());

    // Assign in-bound values to the field and make sure it passes the within-interval check
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f1_data,num_reals,engine,pos_pdf);
    f1.sync_to_dev();
    REQUIRE(interval_check->check(f1));

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f1_data[i] *= -1;
    }
    f1.sync_to_dev();
    REQUIRE(not interval_check->check(f1));
    interval_check->repair(f1);
    REQUIRE(interval_check->check(f1));
  }

  // Check that the values of a field are above a lower bound
  SECTION ("field_lower_bound_check") {
    Field f1(fid);
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    auto lower_bound_check = std::make_shared<FieldLowerBoundCheck>(-1.0);
    REQUIRE(lower_bound_check->can_repair());

    // Assign in-bound values to the field and make sure it passes the lower_bound check
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f1_data[i] = std::numeric_limits<Real>::max() - i*1.0; 
    }
    f1.sync_to_dev();
    REQUIRE(lower_bound_check->check(f1));

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f1_data[i] = -2.0*(i+1);
    }
    f1.sync_to_dev();
    REQUIRE(not lower_bound_check->check(f1));
    lower_bound_check->repair(f1);
    REQUIRE(lower_bound_check->check(f1));
    // Should have repaired to the lower bound:
    f1.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f1_data[i] == -1.0);
    }
  }

  // Check that the values of a field are above below an upper bound
  SECTION ("field_upper_bound_check") {
    Field f1(fid);
    auto upper_bound_check = std::make_shared<FieldUpperBoundCheck>(1.0);
    REQUIRE(upper_bound_check->can_repair());
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    // Assign in-bound values to the field and make sure it passes the upper_bound check
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f1_data[i] = -std::numeric_limits<Real>::max() + i*1.0; 
    }
    f1.sync_to_dev();
    REQUIRE(upper_bound_check->check(f1));

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f1_data[i] = 2.0*(i+1);
    }
    f1.sync_to_dev();
    REQUIRE(not upper_bound_check->check(f1));
    upper_bound_check->repair(f1);
    REQUIRE(upper_bound_check->check(f1));
    // Should have repaired to the upper bound:
    f1.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f1_data[i] == 1.0);
    }
  }

  SECTION ("check_and_repair_wrapper") {
    Field f1(fid);
    f1.allocate_view();
    const int num_reals = f1.get_header().get_alloc_properties().get_num_scalars();

    // Two separate FPC for check and for repair
    constexpr Real ub_check  = 1.0;
    constexpr Real ub_repair = 0.0;
    auto check  = std::make_shared<FieldUpperBoundCheck>(ub_check);
    auto repair = std::make_shared<FieldUpperBoundCheck>(ub_repair);

    auto check_and_repair = std::make_shared<CheckAndRepairWrapper>(check, repair);
    REQUIRE(check_and_repair->can_repair());

    // Assign out-of-bound values to the field, and ensure check fails
    auto f1_data = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f1_data[i] = 2.0;
    }
    f1.sync_to_dev();
    REQUIRE(not check_and_repair->check(f1));

    // Repair the field, and make sure the field values match ub_repair
    check_and_repair->repair(f1);
    f1.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f1_data[i] == ub_repair);
    }
  }
}

} // anonymous namespace
