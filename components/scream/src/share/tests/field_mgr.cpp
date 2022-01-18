#include <catch2/catch.hpp>

#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/field/field_manager.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/user_provided_grids_manager.hpp"

#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"

namespace {

TEST_CASE("field_mgr", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;
  using FR  = FieldRequest;
  using SL  = std::list<std::string>;
  using Pack = ekat::Pack<Real,8>;

  const int ncols = 4;
  const int nlevs = 7;
  const int subview_dim = 1;
  const int subview_slice = 0;

  std::vector<FieldTag> tags1 = {COL,LEV};
  std::vector<FieldTag> tags2 = {EL,GP,GP};
  std::vector<FieldTag> tags3 = {COL,CMP,LEV};

  std::vector<int> dims1 = {ncols,nlevs};
  std::vector<int> dims2 = {2,4,4};
  std::vector<int> dims3 = {ncols,nlevs+1};
  std::vector<int> dims4 = {ncols,10,nlevs};

  const auto km = 1000*m;

  FID fid1("field_1", {tags1, dims1},  m/s, "phys");
  FID fid2("field_2", {tags1, dims1},  m/s, "phys");
  FID fid3("field_3", {tags1, dims1},  m/s, "phys");
  FID fid4("field_4", {tags3, dims4},  m/s, "phys");
  FID fid5("field_5", {tags1, dims1},  m/s, "phys");

  FID bad1("field_1", {tags2, dims2},  m/s, "dyn");  // Bad grid
  FID bad2("field_2", {tags1, dims1}, km/s, "phys"); // Bad units
  FID bad3("field_2", {tags1, dims3},  m/s, "phys"); // Bad layout

  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid("phys",ncols*comm.size(),nlevs,comm);
  FieldManager<Real> field_mgr(pg);

  // Should not be able to register fields yet
  REQUIRE_THROWS(field_mgr.register_field(FR(fid1)));

  field_mgr.registration_begins();

  // === Valid registration calls === //
  field_mgr.register_field(FR(fid1,"group_1",Pack::n));
  field_mgr.register_field(FR{fid2,"group_2",16});
  field_mgr.register_field(FR{fid3,"group_4"});
  field_mgr.register_field(FR{fid3,SL{"group_1","group_2","group_3"}});
  field_mgr.register_field(FR{fid2,"group_4"});
  field_mgr.register_field(FR{fid4});

  // === Invalid registration calls === //
  REQUIRE_THROWS(field_mgr.register_field(FR{bad1}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));

  field_mgr.registration_ends();

  // Should not be able to register fields anymore
  REQUIRE_THROWS(field_mgr.register_field(FR{fid1}));

  // Check registration is indeed closed
  REQUIRE (field_mgr.repository_state()==RepoState::Closed);
  REQUIRE (field_mgr.size()==5);

  // Get all fields
  auto f1 = field_mgr.get_field(fid1.name());
  auto f2 = field_mgr.get_field(fid2.name());
  auto f3 = field_mgr.get_field(fid3.name());
  auto f4 = field_mgr.get_field(fid4.name());
  auto f5 = field_mgr.get_field(fid5.name());
  REQUIRE_THROWS(field_mgr.get_field("bad")); // Not in the field_mgr
  REQUIRE(f1.get_header().get_identifier()==fid1);

  // Check that the groups names are in the header. While at it, make sure that case insensitive works fine.
  auto has_group = [](const ekat::WeakPtrSet<const FieldGroupInfo>& groups,
                      const std::string& name)->bool {
    for (auto it : groups) {
      if (it.lock()->m_group_name==name) {
        return true;
      }
    }
    return false;
  };
  REQUIRE (has_group(f1.get_header().get_tracking().get_groups_info(),"gRouP_1"));
  REQUIRE (has_group(f2.get_header().get_tracking().get_groups_info(),"Group_2"));
  REQUIRE (has_group(f2.get_header().get_tracking().get_groups_info(),"Group_4"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_1"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_2"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_3"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_4"));

  // Check that the groups in the field_mgr contain the correct fields
  REQUIRE (field_mgr.get_groups_info().count("GROUP_1")==1);
  REQUIRE (field_mgr.get_groups_info().count("GRoup_2")==1);
  REQUIRE (field_mgr.get_groups_info().count("group_3")==1);
  REQUIRE (field_mgr.get_groups_info().count("groUP_4")==1);
  REQUIRE (field_mgr.get_groups_info().count("group_5")==0);
  REQUIRE (field_mgr.get_groups_info().at("group_2")->m_fields_names.size()==2);

  auto g1 = field_mgr.get_groups_info().at("group_1");
  auto g2 = field_mgr.get_groups_info().at("group_2");
  auto g3 = field_mgr.get_groups_info().at("group_3");
  auto g4 = field_mgr.get_groups_info().at("group_4");
  REQUIRE (ekat::contains(g1->m_fields_names,"Field_1"));
  REQUIRE (ekat::contains(g1->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(g2->m_fields_names,"Field_2"));
  REQUIRE (ekat::contains(g2->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(g3->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(g4->m_fields_names,"Field_2"));
  REQUIRE (ekat::contains(g4->m_fields_names,"Field_3"));

  // Check alloc props for f1 and f2 (which requested pack size > 1)
  auto f1_padding = f1.get_header().get_alloc_properties().get_padding();
  auto f2_padding = f2.get_header().get_alloc_properties().get_padding();

  REQUIRE (f1_padding==ekat::PackInfo<Pack::n>::padding(nlevs));
  REQUIRE (f2_padding==ekat::PackInfo<16>::padding(nlevs));

  // Verify f5 is a subfield of f4
  auto f5_ap = f5.get_header().get_alloc_properties();
  REQUIRE (f5_ap.is_subfield());
  REQUIRE (f5_ap.is_dynamic_subfield());
  REQUIRE (f5_ap.get_subview_info().dim_idx==subview_dim);
  REQUIRE (f5_ap.get_subview_info().slice_idx==subview_slice);

  // Fill f4 with random numbers, and verify corresponding subview of f5 gets same values.
  auto engine = setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0,1.0);
  randomize(f5,engine,pdf);
  REQUIRE (views_are_equal(f5,f4.get_component(subview_slice)));
}

TEST_CASE("tracers_bundle", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR  = FieldRequest;

  const int ncols = 4;
  const int nlevs = 7;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {ncols,nlevs};

  const auto nondim = Units::nondimensional();

  const std::string grid_name = "physics";
  FieldIdentifier qv_id("qv", {tags, dims}, nondim, grid_name);
  FieldIdentifier qc_id("qc", {tags, dims}, nondim, grid_name);
  FieldIdentifier qr_id("qr", {tags, dims}, nondim, grid_name);

  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid(grid_name,ncols*comm.size(),nlevs,comm);

  FieldManager<Real> field_mgr(pg);
  field_mgr.registration_begins();
  field_mgr.register_field(FR{qv_id,"tracers"});
  field_mgr.register_field(FR{qc_id,"tracers"});
  field_mgr.register_field(FR{qr_id,"tracers"});
  field_mgr.register_group(GroupRequest("tracers",grid_name,Bundling::Required));
  field_mgr.registration_ends();

  auto qv = field_mgr.get_field(qv_id.name());
  auto qc = field_mgr.get_field(qc_id.name());
  auto qr = field_mgr.get_field(qr_id.name());

  // The field_mgr should have allocated the group bundled
  auto group = field_mgr.get_field_group("tracers");
  REQUIRE (group.m_info->m_bundled);

  const auto& Q_name = group.m_bundle->get_header().get_identifier().name();
  auto Q = field_mgr.get_field(Q_name);

  // The bundled field in the group should match the field we get from the field_mgr
  REQUIRE (Q.equivalent(*group.m_bundle));

  // Check that Q is set as parent for all q's.
  auto qvp = qv.get_header().get_parent().lock();
  auto qcp = qc.get_header().get_parent().lock();
  auto qrp = qr.get_header().get_parent().lock();
  REQUIRE ((qvp!=nullptr && qvp.get()==&Q.get_header()));
  REQUIRE ((qcp!=nullptr && qvp.get()==&Q.get_header()));
  REQUIRE ((qrp!=nullptr && qvp.get()==&Q.get_header()));

  // The indices used for each q to subview Q
  int idx_v, idx_c, idx_r;

  // The idx must be stored
  REQUIRE_NOTHROW (idx_v = group.m_info->m_subview_idx.at("qv"));
  REQUIRE_NOTHROW (idx_c = group.m_info->m_subview_idx.at("qc"));
  REQUIRE_NOTHROW (idx_r = group.m_info->m_subview_idx.at("qr"));

  // All idx must be in [0,2] and must be different
  REQUIRE ((idx_v>=0 && idx_v<3 &&
            idx_c>=0 && idx_c<3 &&
            idx_r>=0 && idx_r<3));
  REQUIRE ((idx_v!=idx_c && idx_v!=idx_r && idx_c!=idx_r));

  // Now fill Q with random values
  auto engine = setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0,1.0);

  randomize(Q,engine,pdf);

  // Check that the same values are in all q's
  Q.sync_to_host();
  auto Qh = Q.get_view<Real***,Host>();
  auto qvh = qv.get_view<Real**,Host>();
  auto qch = qc.get_view<Real**,Host>();
  auto qrh = qr.get_view<Real**,Host>();

  for (int icol=0; icol<ncols; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      REQUIRE (Qh(icol,idx_v,ilev)==qvh(icol,ilev));
      REQUIRE (Qh(icol,idx_c,ilev)==qch(icol,ilev));
      REQUIRE (Qh(icol,idx_r,ilev)==qrh(icol,ilev));
    }
  }

  // Check that the field ptrs stored in the group are the same as the q
  auto qv_ptr = group.m_fields.at("qv");
  auto qc_ptr = group.m_fields.at("qc");
  auto qr_ptr = group.m_fields.at("qr");

  REQUIRE (qv_ptr->equivalent(qv));
  REQUIRE (qc_ptr->equivalent(qc));
  REQUIRE (qr_ptr->equivalent(qr));
}

TEST_CASE("multiple_bundles") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using SL = std::list<std::string>;

  const int ncols = 4;
  const int nlevs = 7;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {ncols,nlevs};

  const auto nondim = Units::nondimensional();

  const std::string grid_name = "physics";
  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid(grid_name,ncols*comm.size(),nlevs,comm);

  FieldIdentifier a_id("a", {tags, dims}, nondim, grid_name);
  FieldIdentifier b_id("b", {tags, dims}, nondim, grid_name);
  FieldIdentifier c_id("c", {tags, dims}, nondim, grid_name);
  FieldIdentifier d_id("d", {tags, dims}, nondim, grid_name);
  FieldIdentifier e_id("e", {tags, dims}, nondim, grid_name);
  FieldIdentifier f_id("f", {tags, dims}, nondim, grid_name);

  FieldRequest a_req(a_id,SL{"group1","group3"});
  FieldRequest b_req(b_id,SL{"group1"});
  FieldRequest c_req(c_id,SL{"group1","group2"});
  FieldRequest d_req(d_id,SL{"group1","group3"});
  FieldRequest e_req(e_id,SL{"group1","group2"});
  FieldRequest f_req(f_id,SL{"group1"});

  GroupRequest g1_req ("group1",grid_name,Bundling::Required);
  GroupRequest g2_req ("group2",grid_name,Bundling::Required);
  // Include all group2 in group3
  GroupRequest g3_req ("group3",grid_name,4,Bundling::Required,DerivationType::Superset,g2_req.name,g2_req.grid);
  // Create group4 as a copy of group2
  GroupRequest g4_req ("group4",grid_name,4,Bundling::Required,DerivationType::Copy,g2_req.name,g2_req.grid);
  // Extend group5 by adding all fields in group1 *except* 'c' and 'd'.
  GroupRequest g5_req ("group5",grid_name,4,Bundling::Preferred,DerivationType::Subset,g1_req.name,g1_req.grid,SL{"c","d"});

  // The above group specs should give the following groups:
  // g1: [a,b,c,d,e,f]
  // g2: [c,e]
  // g3: [a,c,d,e]
  // g4: [c,e]
  // g5: [a,b,e,f]
  // The bundling requests can be accommodated for g1,g2,g3,g4, but not g5.
  // But g5 request is only 'Preferred', so the FM won't error out.
  // The order of fields in the 'encompassing' group is {[c,e],[a,d],[b,f]},
  // where [f1,..,fn] means that the order of those two fields can be anything.
  // The 'block'-reverse of that list is also possible: {[b,f],[a,d],[c,e]}

  FieldManager<Real> field_mgr(pg);
  field_mgr.registration_begins();

  // Register single fields
  field_mgr.register_field(a_req);
  field_mgr.register_field(b_req);
  field_mgr.register_field(c_req);
  field_mgr.register_field(d_req);
  field_mgr.register_field(e_req);
  field_mgr.register_field(f_req);

  // Register groups
  field_mgr.register_group(g1_req);
  field_mgr.register_group(g2_req);
  field_mgr.register_group(g3_req);
  field_mgr.register_group(g4_req);
  field_mgr.register_group(g5_req);

  field_mgr.registration_ends();

  auto g1 = field_mgr.get_field_group(g1_req.name);
  auto g2 = field_mgr.get_field_group(g2_req.name);
  auto g3 = field_mgr.get_field_group(g3_req.name);
  auto g4 = field_mgr.get_field_group(g4_req.name);
  auto g5 = field_mgr.get_field_group(g5_req.name);

  // First 4 groups should be bundled
  REQUIRE (g1.m_info->m_bundled);
  REQUIRE (g2.m_info->m_bundled);
  REQUIRE (g3.m_info->m_bundled);
  REQUIRE (g4.m_info->m_bundled);
  REQUIRE (not g5.m_info->m_bundled);

  // Check that the order of fields in g1 is the expected one
  const auto& fnames = g1.m_info->m_fields_names;
  const auto& f1 = *std::next(fnames.begin(),0);
  const auto& f2 = *std::next(fnames.begin(),1);
  const auto& f3 = *std::next(fnames.begin(),2);
  const auto& f4 = *std::next(fnames.begin(),3);
  const auto& f5 = *std::next(fnames.begin(),4);
  const auto& f6 = *std::next(fnames.begin(),5);
  if (f1=="b" || f1=="f") {
    REQUIRE ( ((f1=="b" && f2=="f") || (f1=="f" && f2=="b")) );
    REQUIRE ( ((f3=="a" && f4=="d") || (f3=="d" && f4=="a")) );
    REQUIRE ( ((f5=="c" && f6=="e") || (f5=="e" && f6=="c")) );
  } else {
    REQUIRE ( ((f1=="c" && f2=="e") || (f1=="e" && f2=="c")) );
    REQUIRE ( ((f3=="a" && f4=="d") || (f3=="d" && f4=="a")) );
    REQUIRE ( ((f5=="b" && f6=="f") || (f5=="f" && f6=="b")) );
  }
}

} // anonymous namespace
