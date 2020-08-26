#include "physics/haero/haero_inputs_initializer.hpp"

namespace scream {

void HaeroInputsInitializer::add_field(const field_type& f) {
  const auto& id = f.get_header().get_identifier();

  m_fields.emplace(id.name(), f);
  m_fields_id.insert(id);
}

void HaeroInputsInitializer::initialize_fields() {
  // Safety check: if we're asked to init anything at all,
  // then we should have been asked to init TODO: XX fields.
  int count = 0;
  // count += m_fields.count("q");
  if (count == 0) {
    return;
  }

  //TODO: EKAT_REQUIRE_MSG(count == XX, "MSG");

  // get device views


  // create host mirrors

  // initialize data

  // deep copy to device

  /* if we are in charge of init-ing <SHARED_FIELD>, init to 0
  if (m_fields.count("<SHARED_FIELD>") == 1) {
    auto fshr = m_fields.at("<SHARED_FIELD>").get_view();
    Kokkos::deep_copy(fshr, Real(0));
  }*/
}

} // namespace scream
