"""Tests for automatic CHEMPOTPY surface discovery."""

import json
import unittest
from pathlib import Path

from chempotpy.core.registry import (
    discover_surface_sources,
    duplicate_surface_sources,
)


class RegistryTests( unittest.TestCase ):
    def test_discovers_library_sources( self ) -> None:
        sources = discover_surface_sources()
        self.assertGreater(len(sources), 300)
        self.assertGreater(len({source.system_name for source in sources}), 80)
        self.assertTrue(all(source.path.is_file() for source in sources))
        fixture = json.loads(
            (
                Path(__file__).resolve().parent
                / "data"
                / "all_surface_references.json"
            ).read_text()
        )["surfaces"]
        self.assertEqual(
            sum(source.entry_point == "pes" for source in sources),
            sum(record["entry_point"] == "pes" for record in fixture.values()),
        )
        self.assertEqual(
            sum(source.entry_point == "dpem" for source in sources),
            sum(record["entry_point"] == "dpem" for record in fixture.values()),
        )
        support_sources = {
            source.path.name for source in sources if not source.is_surface
        }
        self.assertEqual(
            support_sources,
            {
                "ch4_t8.f90",
                "chempotpy_arno_spline_parameters.f90",
                "chempotpy_c2h2n2o_pipnn_parameters.f90",
                "chempotpy_c2h2o_pipnn_parameters.f90",
                "chempotpy_c2h4o4_invariants.f90",
                "chempotpy_c2h4o4_parameters.f90",
                "chempotpy_c3h7no_pip_parameters.f90",
                "chempotpy_c2n_pipnn_ap_parameters.f90",
                "chempotpy_c2n_pipnn_app_parameters.f90",
                "chempotpy_co2_1ap_parameters.f90",
                "chempotpy_co2_1app_parameters.f90",
                "chempotpy_co2_21ap_parameters.f90",
                "chempotpy_co2_3ap_parameters.f90",
                "chempotpy_co2_3app_parameters.f90",
                "chempotpy_co2_asymptotic_parameters.f90",
                "chempotpy_co2_rkhs_runtime.f90",
                "chempotpy_ch2oh_2state_parameters.f90",
                "chempotpy_ch2oh_3state_parameters.f90",
                "chempotpy_ch2oh_surfgen_basis.f90",
                "chempotpy_ch2oh_surfgen_runtime.f90",
                "chempotpy_ch2o2_pipnn_parameters.f90",
                "chempotpy_ch2o2_local_pipnn_parameters.f90",
                "chempotpy_ch2o_pipnn_parameters.f90",
                "chempotpy_ch4_t8_parameters.f90",
                "chempotpy_ch4ni_ni111_parameters.f90",
                "chempotpy_ch4ni_runtime.f90",
                "chempotpy_ch4ocl_pipnn_parameters.f90",
                "chempotpy_coco_pipnn_parameters.f90",
                "chempotpy_ch4of_pipnn_parameters.f90",
                "chempotpy_clh2o_pipnn_parameters.f90",
                "chempotpy_f2h2_jbkkl_parameters.f90",
                "chempotpy_fh2o_pipnn_parameters.f90",
                "chempotpy_h2co_pip_parameters.f90",
                "chempotpy_n4_pipnn_parameters.f90",
                "chempotpy_hobr_3app_parameters.f90",
                "chempotpy_h3_2002_basis.f90",
                "chempotpy_h3cl_pipnn_parameters.f90",
                "chempotpy_h3clo_pipnn_parameters.f90",
                "chempotpy_hcl_h2o_pip_basis.f90",
                "chempotpy_hcl_h2o_pip_parameters.f90",
                "chempotpy_h3s_pipnn_parameters.f90",
                "chempotpy_h4o2_pip_parameters.f90",
                "chempotpy_h7p_pip_basis.f90",
                "chempotpy_h7p_pip_parameters.f90",
                "chempotpy_hni_eam5_parameters.f90",
                "chempotpy_hni_eam67_parameters.f90",
                "chempotpy_hni_eam67_runtime.f90",
                "chempotpy_hni_ni111_parameters.f90",
                "chempotpy_hcloh_pipnn_parameters.f90",
                "chempotpy_ho3_pipnn_parameters.f90",
                "chempotpy_k2rb2_pipnn_parameters.f90",
                "chempotpy_malonaldehyde_mht_parameters.f90",
                "chempotpy_n2hocp_pipnn_parameters.f90",
                "chempotpy_n4_bdtc_parameters.f90",
                "chempotpy_nh3_surfgen_parameters.f90",
                "chempotpy_nh3_surfgen_runtime.f90",
                "chempotpy_no2_2ap_pipnn_parameters.f90",
                "chempotpy_no2_4ap_pipnn_parameters.f90",
                "chempotpy_no2_6ap_pipnn_parameters.f90",
                "chempotpy_o2h3_pipnn_parameters.f90",
            },
        )

    def test_reports_competing_fixed_and_free_form_sources( self ) -> None:
        duplicates = duplicate_surface_sources(discover_surface_sources())
        self.assertEqual(duplicates, {})


if __name__ == "__main__":
    unittest.main()
