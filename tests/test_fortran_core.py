"""Compile and exercise the bundled modern Fortran core modules."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class FortranCoreTests( unittest.TestCase ):
    def test_core_and_migrated_surface_support_quad_precision_when_available(
        self,
    ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core"
        with tempfile.TemporaryDirectory(
            prefix="chempotpy-quad-precision-"
        ) as tmp:
            directory = Path(tmp)
            probe = directory / "probe.f90"
            probe_executable = directory / "probe"
            probe.write_text("""\
program probe_quad_kind
    implicit none
    print *, selected_real_kind(33, 4931)
end program probe_quad_kind
""")
            subprocess.run(
                [compiler, str(probe), "-o", str(probe_executable)],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            quad_kind = int(subprocess.run(
                [str(probe_executable)],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            ).stdout.strip())
            if quad_kind <= 0:
                self.skipTest("compiler has no requested quadruple real kind")
            kinds = directory / "chempotpy_kinds.f90"
            kinds.write_text(
                (core / "chempotpy_kinds.f90").read_text().replace(
                    "wp = dp", "wp = qp"
                )
            )
            driver = directory / "driver.f90"
            executable = directory / "test_quad_precision"
            driver.write_text("""\
program test_quad_precision
    use chempotpy_aln_er_pes, only: evaluate
    use chempotpy_kinds, only: wp
    use chempotpy_units, only: angstrom_to_bohr, bohr_to_angstrom
    implicit none
    real(wp) :: coordinates(10000,3), energy(1)
    real(wp) :: gradient(1,2,3), coupling(1,1,2,3)
    real(wp) :: round_trip
    if( precision(0.0_wp) < 33 .or. range(0.0_wp) < 4931 ) &
        error stop 'wp did not select requested quadruple precision'
    round_trip = angstrom_to_bohr*bohr_to_angstrom
    if( abs(round_trip-1.0_wp) > epsilon(1.0_wp) ) &
        error stop 'quadruple-precision unit round trip failed'
    coordinates = 0.0_wp
    coordinates(2,1) = 2.8_wp
    call evaluate(coordinates, 2, 0, energy, gradient, coupling)
    if( energy(1) /= energy(1) ) &
        error stop 'quadruple-precision PES evaluation returned NaN'
end program test_quad_precision
""")
            subprocess.run(
                [
                    compiler,
                    str(kinds),
                    str(core / "chempotpy_physical_constants.f90"),
                    str(core / "chempotpy_units.f90"),
                    str(core / "chempotpy_geometry.f90"),
                    str(repository / "chempotpy" / "Aln" / "Aln_ER.f90"),
                    str(driver),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            subprocess.run([str(executable)], check=True, cwd=directory)

    def test_units_geometry_and_status_contract( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core"
        source_names = (
            "chempotpy_kinds.f90",
            "chempotpy_physical_constants.f90",
            "chempotpy_units.f90",
            "chempotpy_legacy_units.f90",
            "chempotpy_geometry.f90",
            "chempotpy_internal_coordinates.f90",
            "chempotpy_legacy_coordinates.f90",
            "chempotpy_status.f90",
        )
        driver_source = """
program test_fortran_core
    use chempotpy_geometry, only: all_pair_distance_jacobian, all_pair_distances, &
        all_pair_distances_and_jacobian, cartesian_distance, cosine_between_vectors, dot_product_3d, fill_distance_matrix, &
        flatten_cartesian_coordinates, ordered_pair_distance_jacobian, transform_internal_gradient_to_cartesian, &
        triatomic_distance_jacobian
    use chempotpy_internal_coordinates, only: finalize_legacy_triatomic_surface, &
        prepare_legacy_triatomic_coordinates
    use chempotpy_kinds, only: dp, real32, real64, sp, wp
    use chempotpy_physical_constants, only: codata_adjustment_year
    use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates, &
        legacy_triatomic_distances_to_bohr, triatomic_distances_to_bohr
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, &
        legacy_coordinates_to_bohr, legacy_energy_to_electronvolt, &
        legacy_gradient_to_electronvolt_per_angstrom, &
        legacy_hartree_to_electronvolt
    use chempotpy_status, only: status_derivative_unavailable, &
        status_invalid_input, status_success
    use chempotpy_units, only: angstrom_to_bohr, bohr_to_angstrom, &
        electronvolt_to_hartree, electronvolt_to_inverse_centimeter, &
        hartree_per_bohr_to_ev_per_angstrom, hartree_to_electronvolt, &
        hartree_to_inverse_centimeter, hartree_to_kilocalorie_per_mole
    implicit none

    real(wp) :: first_point(3), second_point(3), coordinates(3,3)
    real(wp) :: distances(3,3), converted(1,3), flat_coordinates(9)
    real(wp) :: coordinate_jacobian(2,3), internal_gradient(2)
    real(wp) :: internal_gradient_matrix(2,2), cartesian_gradient_matrix(2,3)
    real(wp) :: internal_gradient_tensor(2,2,2), cartesian_gradient_tensor(2,2,3)
    real(wp) :: cartesian_gradient(3)
    real(wp) :: potlib_coordinates(9)
    real(wp) :: pair_distances(3)
    real(wp) :: triatomic_jacobian(3,9)
    real(wp) :: four_atom_coordinates(4,3), six_distances(6)
    real(wp) :: six_distance_jacobian(6,12), flat_four_atom_coordinates(12)
    real(wp) :: surface_energy(1), surface_gradient(1,3,3)
    real(wp) :: surface_coupling(1,1,3,3)
    real(wp) :: surface_energy_hartree, surface_internal_gradient(3)
    integer      :: conversion_status

    if( sp /= real32 ) error stop 'single-precision kind is not binary32'
    if( dp /= real64 ) error stop 'double-precision kind is not binary64'
    if( wp /= dp ) error stop 'default working precision is not double precision'
    if( precision(0.0_wp) < 15 .or. range(0.0_wp) < 307 ) &
        error stop 'working precision contract is too small'
    if( codata_adjustment_year /= 2022 ) error stop 'CODATA version changed'
    if( bohr_to_angstrom /= 0.529177210544_wp ) &
        error stop 'canonical length conversion changed'
    if( hartree_to_electronvolt /= 27.211386245981_wp ) &
        error stop 'canonical energy conversion changed'
    if( abs(hartree_per_bohr_to_ev_per_angstrom-51.42206751119799_wp) &
        > 4.0_wp*epsilon(hartree_per_bohr_to_ev_per_angstrom)* &
        abs(hartree_per_bohr_to_ev_per_angstrom) ) &
        error stop 'canonical gradient conversion changed'
    if( abs(hartree_to_kilocalorie_per_mole-627.5094740628975_wp) &
        > 4.0_wp*epsilon(hartree_to_kilocalorie_per_mole)* &
        abs(hartree_to_kilocalorie_per_mole) ) &
        error stop 'canonical molar-energy conversion changed'
    if( abs(hartree_to_inverse_centimeter-219474.63136314_wp) &
        > 4.0_wp*epsilon(hartree_to_inverse_centimeter)* &
        abs(hartree_to_inverse_centimeter) ) &
        error stop 'canonical Hartree wavenumber changed'
    if( abs(electronvolt_to_inverse_centimeter-8065.543937349212_wp) &
        > 4.0_wp*epsilon(electronvolt_to_inverse_centimeter)* &
        abs(electronvolt_to_inverse_centimeter) ) &
        error stop 'canonical electronvolt wavenumber changed'
    if( abs(angstrom_to_bohr*bohr_to_angstrom-1.0_wp) &
        > epsilon(1.0_wp) ) error stop 'length conversion is inconsistent'
    if( abs(electronvolt_to_hartree*hartree_to_electronvolt-1.0_wp) &
        > epsilon(1.0_wp) ) error stop 'energy conversion is inconsistent'

    if( legacy_bohr_to_angstrom /= 0.52917718887329101562_wp ) &
        error stop 'legacy length rounding changed'
    if( legacy_hartree_to_electronvolt /= &
        27.211385726928710938_wp ) &
        error stop 'legacy energy rounding changed'
    call legacy_coordinates_to_bohr( &
        reshape([legacy_bohr_to_angstrom, 0.0_wp, 0.0_wp], [1,3]), &
        converted )
    if( converted(1,1) /= 1.0_wp ) &
        error stop 'legacy coordinate conversion failed'
    if( legacy_energy_to_electronvolt(1.0_wp) /= &
        legacy_hartree_to_electronvolt ) &
        error stop 'legacy energy conversion failed'

    first_point = [0.0_wp, 0.0_wp, 0.0_wp]
    second_point = [3.0_wp, 4.0_wp, 12.0_wp]
    if( cartesian_distance(first_point,second_point) /= 13.0_wp ) &
        error stop 'Cartesian distance failed'
    if( dot_product_3d([1.0_wp,2.0_wp,3.0_wp], &
        [4.0_wp,5.0_wp,6.0_wp]) /= 32.0_wp ) &
        error stop 'dot product failed'
    if( abs(cosine_between_vectors([1.0_wp,0.0_wp,0.0_wp], &
        [0.0_wp,1.0_wp,0.0_wp])) > epsilon(1.0_wp) ) &
        error stop 'angle cosine failed'
    coordinates = 0.0_wp
    coordinates(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    coordinates(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    coordinates(3,:) = [7.0_wp, 8.0_wp, 9.0_wp]
    call flatten_cartesian_coordinates(coordinates, flat_coordinates)
    if( any(flat_coordinates /= [1.0_wp, 2.0_wp, 3.0_wp, &
        4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, &
        9.0_wp]) ) error stop 'Cartesian flattening order changed'
    call cartesian_to_potlib_coordinates( &
        coordinates, 3, 1, potlib_coordinates, conversion_status )
    if( conversion_status /= 0 .or. any(potlib_coordinates /= flat_coordinates) ) &
        error stop 'POTLIB Cartesian flattening changed'
    coordinate_jacobian(1,:) = [1.0_wp, 2.0_wp, 3.0_wp]
    coordinate_jacobian(2,:) = [4.0_wp, 5.0_wp, 6.0_wp]
    internal_gradient = [2.0_wp, 3.0_wp]
    call transform_internal_gradient_to_cartesian( &
        internal_gradient, coordinate_jacobian, cartesian_gradient )
    if( any(cartesian_gradient /= [14.0_wp, 19.0_wp, 24.0_wp]) ) &
        error stop 'vector chain-rule transformation failed'
    internal_gradient_matrix(1,:) = internal_gradient
    internal_gradient_matrix(2,:) = [-1.0_wp, 4.0_wp]
    call transform_internal_gradient_to_cartesian( &
        internal_gradient_matrix, coordinate_jacobian, cartesian_gradient_matrix )
    if( any(cartesian_gradient_matrix(1,:) /= cartesian_gradient) .or. &
        any(cartesian_gradient_matrix(2,:) /= [15.0_wp, 18.0_wp, 21.0_wp]) ) &
        error stop 'matrix chain-rule transformation failed'
    internal_gradient_tensor = 0.0_wp
    internal_gradient_tensor(1,2,:) = internal_gradient
    call transform_internal_gradient_to_cartesian( &
        internal_gradient_tensor, coordinate_jacobian, cartesian_gradient_tensor )
    if( any(cartesian_gradient_tensor(1,2,:) /= cartesian_gradient) ) &
        error stop 'tensor chain-rule transformation failed'
    coordinates = 0.0_wp
    coordinates(2,1) = 3.0_wp
    coordinates(3,2) = 4.0_wp
    call fill_distance_matrix( coordinates, distances )
    if( distances(1,2) /= 3.0_wp &
        .or. distances(1,3) /= 4.0_wp &
        .or. distances(2,3) /= 5.0_wp ) &
        error stop 'distance matrix failed'
    if( any(distances /= transpose(distances)) ) &
        error stop 'distance matrix is not symmetric'
    call fill_distance_matrix(coordinates(:,1), coordinates(:,2), coordinates(:,3), distances)
    if( distances(1,2) /= 3.0_wp .or. distances(1,3) /= 4.0_wp .or. &
        distances(2,3) /= 5.0_wp ) error stop 'component-array distance matrix failed'
    coordinates = 0.0_wp
    coordinates(2,1) = 3.0_wp*legacy_bohr_to_angstrom
    coordinates(3,1) = 3.0_wp*legacy_bohr_to_angstrom
    coordinates(3,2) = 4.0_wp*legacy_bohr_to_angstrom
    call legacy_triatomic_distances_to_bohr(coordinates, pair_distances)
    if( any(abs(pair_distances-[3.0_wp, 4.0_wp, 5.0_wp]) &
        > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'legacy triatomic distances failed'
    coordinates(2,1) = 3.0_wp*bohr_to_angstrom
    coordinates(3,1) = 3.0_wp*bohr_to_angstrom
    coordinates(3,2) = 4.0_wp*bohr_to_angstrom
    call triatomic_distances_to_bohr(coordinates, pair_distances)
    if( any(abs(pair_distances-[3.0_wp, 4.0_wp, 5.0_wp]) &
        > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'canonical triatomic distances failed'
    call triatomic_distance_jacobian(coordinates, pair_distances, triatomic_jacobian)
    if( any(abs(pair_distances-bohr_to_angstrom*[3.0_wp, 4.0_wp, 5.0_wp]) &
        > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'triatomic distance Jacobian distances failed'
    if( any(abs(triatomic_jacobian(1,:)-[-1.0_wp, 0.0_wp, 0.0_wp, &
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]) &
        > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'triatomic distance Jacobian orientation failed'
    four_atom_coordinates = 0.0_wp
    four_atom_coordinates(2,1) = 3.0_wp
    four_atom_coordinates(3,2) = 4.0_wp
    four_atom_coordinates(4,3) = 12.0_wp
    call all_pair_distances(four_atom_coordinates, six_distances)
    if( any(six_distances /= [3.0_wp, 4.0_wp, 12.0_wp, 5.0_wp, &
        sqrt(153.0_wp), sqrt(160.0_wp)]) ) error stop 'distance-only matrix operation failed'
    call flatten_cartesian_coordinates(four_atom_coordinates, flat_four_atom_coordinates)
    call all_pair_distances(flat_four_atom_coordinates, six_distances)
    if( any(six_distances /= [3.0_wp, 4.0_wp, 12.0_wp, 5.0_wp, &
        sqrt(153.0_wp), sqrt(160.0_wp)]) ) error stop 'distance-only vector operation failed'
    call all_pair_distances_and_jacobian(four_atom_coordinates, six_distances, six_distance_jacobian)
    if( any(six_distances /= [3.0_wp, 4.0_wp, 12.0_wp, 5.0_wp,  &
        sqrt(153.0_wp), sqrt(160.0_wp)]) ) error stop 'all-pair matrix distances failed'
    call all_pair_distances_and_jacobian(flat_four_atom_coordinates, six_distances, six_distance_jacobian)
    if( any(six_distances /= [3.0_wp, 4.0_wp, 12.0_wp, 5.0_wp, &
        sqrt(153.0_wp), sqrt(160.0_wp)]) ) error stop 'all-pair vector distances failed'
    if( any(abs(six_distance_jacobian(1,:)-[-1.0_wp, 0.0_wp, 0.0_wp, &
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
        0.0_wp, 0.0_wp, 0.0_wp]) > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'all-pair Jacobian orientation failed'
    six_distance_jacobian = 0.0_wp
    call all_pair_distance_jacobian(flat_four_atom_coordinates, six_distances, six_distance_jacobian)
    if( any(abs(six_distance_jacobian(1,:)-[-1.0_wp, 0.0_wp, 0.0_wp, &
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
        0.0_wp, 0.0_wp, 0.0_wp]) > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'supplied-distance all-pair Jacobian failed'
    call ordered_pair_distance_jacobian(flat_four_atom_coordinates, [1, 3, 1], [2, 4, 4], &
        six_distances([1, 6, 3]), six_distance_jacobian(1:3,:))
    if( any(abs(six_distance_jacobian(2,7:12)-[0.0_wp, 4.0_wp/sqrt(160.0_wp), &
        -12.0_wp/sqrt(160.0_wp), 0.0_wp, -4.0_wp/sqrt(160.0_wp), &
        12.0_wp/sqrt(160.0_wp)]) > 8.0_wp*epsilon(1.0_wp)) ) &
        error stop 'ordered-pair Jacobian failed'
    call prepare_legacy_triatomic_coordinates(coordinates, pair_distances, triatomic_jacobian)
    call test_triatomic_kernel(pair_distances, 1, surface_energy_hartree, surface_internal_gradient)
    call finalize_legacy_triatomic_surface( &
        1, surface_energy_hartree, surface_internal_gradient, triatomic_jacobian, &
        surface_energy, surface_gradient, surface_coupling)
    if( surface_energy(1) /= 12.0_wp*hartree_to_electronvolt ) &
        error stop 'internal-coordinate driver energy conversion failed'
    if( any(abs(surface_gradient(1,:,:)-hartree_per_bohr_to_ev_per_angstrom*reshape([ &
        -2.8_wp, 1.0_wp, 1.8_wp, -2.4_wp, -2.0_wp, 4.4_wp, &
        0.0_wp, 0.0_wp, 0.0_wp], [3,3])) > &
        64.0_wp*epsilon(1.0_wp)*hartree_per_bohr_to_ev_per_angstrom*5.0_wp) ) &
        error stop 'internal-coordinate driver gradient transformation failed'
    if( any(surface_coupling /= 0.0_wp) ) &
        error stop 'internal-coordinate driver coupling initialization failed'

    if( status_success /= 0 .or. status_derivative_unavailable /= 1 &
        .or. status_invalid_input /= 2 ) error stop 'status contract changed'
    contains
    subroutine test_triatomic_kernel( distances_bohr, derivative_order, energy_hartree, gradient_hartree_per_bohr )
        real(wp), intent(in)  :: distances_bohr(3)
        integer, intent(in)       :: derivative_order
        real(wp), intent(out) :: energy_hartree, gradient_hartree_per_bohr(3)
        energy_hartree = sum(distances_bohr)
        gradient_hartree_per_bohr = 0.0_wp
        if( derivative_order == 1 ) gradient_hartree_per_bohr = [1.0_wp, 2.0_wp, 3.0_wp]
    end subroutine test_triatomic_kernel
end program test_fortran_core
"""
        with tempfile.TemporaryDirectory(prefix="chempotpy-fortran-core-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_fortran_core"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    *(str(core / name) for name in source_names),
                    str(driver),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    unittest.main()
