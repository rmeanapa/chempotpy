"""Compile and exercise the bundled modern Fortran core modules."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class FortranCoreTests( unittest.TestCase ):
    def test_units_geometry_and_status_contract( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core"
        source_names = (
            "chempotpy_kinds.f90",
            "chempotpy_units.f90",
            "chempotpy_legacy_units.f90",
            "chempotpy_geometry.f90",
            "chempotpy_internal_coordinates.f90",
            "chempotpy_legacy_coordinates.f90",
            "chempotpy_h3_basis.f90",
            "chempotpy_status.f90",
        )
        driver_source = """
program test_fortran_core
    use chempotpy_geometry, only: all_pair_distance_jacobian, all_pair_distances, &
        all_pair_distances_and_jacobian, cartesian_distance, cosine_between_vectors, dot_product_3d, fill_distance_matrix, &
        flatten_cartesian_coordinates, ordered_pair_distance_jacobian, transform_internal_gradient_to_cartesian, &
        triatomic_distance_jacobian
    use chempotpy_h3_basis, only: get_h3_basis_term, primary_term_count, &
        secondary_term_count
    use chempotpy_internal_coordinates, only: finalize_legacy_triatomic_surface, &
        prepare_legacy_triatomic_coordinates
    use chempotpy_kinds, only: real64
    use chempotpy_legacy_coordinates, only: cartesian_to_potlib_coordinates, &
        legacy_triatomic_distances_to_bohr
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, &
        legacy_coordinates_to_bohr, legacy_energy_to_electronvolt, &
        legacy_gradient_to_electronvolt_per_angstrom, &
        legacy_hartree_to_electronvolt
    use chempotpy_status, only: status_derivative_unavailable, &
        status_invalid_input, status_success
    use chempotpy_units, only: angstrom_to_bohr, bohr_to_angstrom, &
        electronvolt_to_hartree, hartree_to_electronvolt
    implicit none

    real(real64) :: first_point(3), second_point(3), coordinates(3,3)
    real(real64) :: distances(3,3), converted(1,3), flat_coordinates(9)
    real(real64) :: coordinate_jacobian(2,3), internal_gradient(2)
    real(real64) :: internal_gradient_matrix(2,2), cartesian_gradient_matrix(2,3)
    real(real64) :: internal_gradient_tensor(2,2,2), cartesian_gradient_tensor(2,2,3)
    real(real64) :: cartesian_gradient(3)
    real(real64) :: potlib_coordinates(9)
    real(real64) :: pair_distances(3)
    real(real64) :: triatomic_jacobian(3,9)
    real(real64) :: four_atom_coordinates(4,3), six_distances(6)
    real(real64) :: six_distance_jacobian(6,12), flat_four_atom_coordinates(12)
    real(real64) :: surface_energy(1), surface_gradient(1,3,3)
    real(real64) :: surface_coupling(1,1,3,3)
    real(real64) :: surface_energy_hartree, surface_internal_gradient(3)
    integer      :: conversion_status
    integer      :: coefficient_index, first_power, second_power, third_power

    if( bohr_to_angstrom /= 0.529177210903_real64 ) &
        error stop 'canonical length conversion changed'
    if( hartree_to_electronvolt /= 27.211386245988_real64 ) &
        error stop 'canonical energy conversion changed'
    if( abs(angstrom_to_bohr*bohr_to_angstrom-1.0_real64) &
        > epsilon(1.0_real64) ) error stop 'length conversion is inconsistent'
    if( abs(electronvolt_to_hartree*hartree_to_electronvolt-1.0_real64) &
        > epsilon(1.0_real64) ) error stop 'energy conversion is inconsistent'

    if( legacy_bohr_to_angstrom /= 0.52917718887329101562_real64 ) &
        error stop 'legacy length rounding changed'
    if( legacy_hartree_to_electronvolt /= &
        27.211385726928710938_real64 ) &
        error stop 'legacy energy rounding changed'
    call legacy_coordinates_to_bohr( &
        reshape([legacy_bohr_to_angstrom, 0.0_real64, 0.0_real64], [1,3]), &
        converted )
    if( converted(1,1) /= 1.0_real64 ) &
        error stop 'legacy coordinate conversion failed'
    if( legacy_energy_to_electronvolt(1.0_real64) /= &
        legacy_hartree_to_electronvolt ) &
        error stop 'legacy energy conversion failed'

    first_point = [0.0_real64, 0.0_real64, 0.0_real64]
    second_point = [3.0_real64, 4.0_real64, 12.0_real64]
    if( cartesian_distance(first_point,second_point) /= 13.0_real64 ) &
        error stop 'Cartesian distance failed'
    if( dot_product_3d([1.0_real64,2.0_real64,3.0_real64], &
        [4.0_real64,5.0_real64,6.0_real64]) /= 32.0_real64 ) &
        error stop 'dot product failed'
    if( abs(cosine_between_vectors([1.0_real64,0.0_real64,0.0_real64], &
        [0.0_real64,1.0_real64,0.0_real64])) > epsilon(1.0_real64) ) &
        error stop 'angle cosine failed'
    coordinates = 0.0_real64
    coordinates(1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
    coordinates(2,:) = [4.0_real64, 5.0_real64, 6.0_real64]
    coordinates(3,:) = [7.0_real64, 8.0_real64, 9.0_real64]
    call flatten_cartesian_coordinates(coordinates, flat_coordinates)
    if( any(flat_coordinates /= [1.0_real64, 2.0_real64, 3.0_real64, &
        4.0_real64, 5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64, &
        9.0_real64]) ) error stop 'Cartesian flattening order changed'
    call cartesian_to_potlib_coordinates( &
        coordinates, 3, 1, potlib_coordinates, conversion_status )
    if( conversion_status /= 0 .or. any(potlib_coordinates /= flat_coordinates) ) &
        error stop 'POTLIB Cartesian flattening changed'
    coordinate_jacobian(1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
    coordinate_jacobian(2,:) = [4.0_real64, 5.0_real64, 6.0_real64]
    internal_gradient = [2.0_real64, 3.0_real64]
    call transform_internal_gradient_to_cartesian( &
        internal_gradient, coordinate_jacobian, cartesian_gradient )
    if( any(cartesian_gradient /= [14.0_real64, 19.0_real64, 24.0_real64]) ) &
        error stop 'vector chain-rule transformation failed'
    internal_gradient_matrix(1,:) = internal_gradient
    internal_gradient_matrix(2,:) = [-1.0_real64, 4.0_real64]
    call transform_internal_gradient_to_cartesian( &
        internal_gradient_matrix, coordinate_jacobian, cartesian_gradient_matrix )
    if( any(cartesian_gradient_matrix(1,:) /= cartesian_gradient) .or. &
        any(cartesian_gradient_matrix(2,:) /= [15.0_real64, 18.0_real64, 21.0_real64]) ) &
        error stop 'matrix chain-rule transformation failed'
    internal_gradient_tensor = 0.0_real64
    internal_gradient_tensor(1,2,:) = internal_gradient
    call transform_internal_gradient_to_cartesian( &
        internal_gradient_tensor, coordinate_jacobian, cartesian_gradient_tensor )
    if( any(cartesian_gradient_tensor(1,2,:) /= cartesian_gradient) ) &
        error stop 'tensor chain-rule transformation failed'
    coordinates = 0.0_real64
    coordinates(2,1) = 3.0_real64
    coordinates(3,2) = 4.0_real64
    call fill_distance_matrix( coordinates, distances )
    if( distances(1,2) /= 3.0_real64 &
        .or. distances(1,3) /= 4.0_real64 &
        .or. distances(2,3) /= 5.0_real64 ) &
        error stop 'distance matrix failed'
    if( any(distances /= transpose(distances)) ) &
        error stop 'distance matrix is not symmetric'
    call fill_distance_matrix(coordinates(:,1), coordinates(:,2), coordinates(:,3), distances)
    if( distances(1,2) /= 3.0_real64 .or. distances(1,3) /= 4.0_real64 .or. &
        distances(2,3) /= 5.0_real64 ) error stop 'component-array distance matrix failed'
    coordinates = 0.0_real64
    coordinates(2,1) = 3.0_real64*legacy_bohr_to_angstrom
    coordinates(3,1) = 3.0_real64*legacy_bohr_to_angstrom
    coordinates(3,2) = 4.0_real64*legacy_bohr_to_angstrom
    call legacy_triatomic_distances_to_bohr(coordinates, pair_distances)
    if( any(abs(pair_distances-[3.0_real64, 4.0_real64, 5.0_real64]) &
        > 8.0_real64*epsilon(1.0_real64)) ) &
        error stop 'legacy triatomic distances failed'
    call triatomic_distance_jacobian(coordinates, pair_distances, triatomic_jacobian)
    if( any(abs(pair_distances-legacy_bohr_to_angstrom*[3.0_real64, 4.0_real64, 5.0_real64]) &
        > 8.0_real64*epsilon(1.0_real64)) ) &
        error stop 'triatomic distance Jacobian distances failed'
    if( any(abs(triatomic_jacobian(1,:)-[-1.0_real64, 0.0_real64, 0.0_real64, &
        1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]) &
        > 8.0_real64*epsilon(1.0_real64)) ) &
        error stop 'triatomic distance Jacobian orientation failed'
    four_atom_coordinates = 0.0_real64
    four_atom_coordinates(2,1) = 3.0_real64
    four_atom_coordinates(3,2) = 4.0_real64
    four_atom_coordinates(4,3) = 12.0_real64
    call all_pair_distances(four_atom_coordinates, six_distances)
    if( any(six_distances /= [3.0_real64, 4.0_real64, 12.0_real64, 5.0_real64, &
        sqrt(153.0_real64), sqrt(160.0_real64)]) ) error stop 'distance-only matrix operation failed'
    call flatten_cartesian_coordinates(four_atom_coordinates, flat_four_atom_coordinates)
    call all_pair_distances(flat_four_atom_coordinates, six_distances)
    if( any(six_distances /= [3.0_real64, 4.0_real64, 12.0_real64, 5.0_real64, &
        sqrt(153.0_real64), sqrt(160.0_real64)]) ) error stop 'distance-only vector operation failed'
    call all_pair_distances_and_jacobian(four_atom_coordinates, six_distances, six_distance_jacobian)
    if( any(six_distances /= [3.0_real64, 4.0_real64, 12.0_real64, 5.0_real64,  &
        sqrt(153.0_real64), sqrt(160.0_real64)]) ) error stop 'all-pair matrix distances failed'
    call all_pair_distances_and_jacobian(flat_four_atom_coordinates, six_distances, six_distance_jacobian)
    if( any(six_distances /= [3.0_real64, 4.0_real64, 12.0_real64, 5.0_real64, &
        sqrt(153.0_real64), sqrt(160.0_real64)]) ) error stop 'all-pair vector distances failed'
    if( any(abs(six_distance_jacobian(1,:)-[-1.0_real64, 0.0_real64, 0.0_real64, &
        1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        0.0_real64, 0.0_real64, 0.0_real64]) > 8.0_real64*epsilon(1.0_real64)) ) &
        error stop 'all-pair Jacobian orientation failed'
    six_distance_jacobian = 0.0_real64
    call all_pair_distance_jacobian(flat_four_atom_coordinates, six_distances, six_distance_jacobian)
    if( any(abs(six_distance_jacobian(1,:)-[-1.0_real64, 0.0_real64, 0.0_real64, &
        1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
        0.0_real64, 0.0_real64, 0.0_real64]) > 8.0_real64*epsilon(1.0_real64)) ) &
        error stop 'supplied-distance all-pair Jacobian failed'
    call ordered_pair_distance_jacobian(flat_four_atom_coordinates, [1, 3, 1], [2, 4, 4], &
        six_distances([1, 6, 3]), six_distance_jacobian(1:3,:))
    if( any(abs(six_distance_jacobian(2,7:12)-[0.0_real64, 4.0_real64/sqrt(160.0_real64), &
        -12.0_real64/sqrt(160.0_real64), 0.0_real64, -4.0_real64/sqrt(160.0_real64), &
        12.0_real64/sqrt(160.0_real64)]) > 8.0_real64*epsilon(1.0_real64)) ) &
        error stop 'ordered-pair Jacobian failed'
    call prepare_legacy_triatomic_coordinates(coordinates, pair_distances, triatomic_jacobian)
    call test_triatomic_kernel(pair_distances, 1, surface_energy_hartree, surface_internal_gradient)
    call finalize_legacy_triatomic_surface( &
        1, surface_energy_hartree, surface_internal_gradient, triatomic_jacobian, &
        surface_energy, surface_gradient, surface_coupling)
    if( surface_energy(1) /= 12.0_real64*legacy_hartree_to_electronvolt ) &
        error stop 'internal-coordinate driver energy conversion failed'
    if( any(abs(surface_gradient(1,:,:)-legacy_gradient_to_electronvolt_per_angstrom*reshape([ &
        -2.8_real64, 1.0_real64, 1.8_real64, -2.4_real64, -2.0_real64, 4.4_real64, &
        0.0_real64, 0.0_real64, 0.0_real64], [3,3])) > 64.0_real64*epsilon(1.0_real64)) ) &
        error stop 'internal-coordinate driver gradient transformation failed'
    if( any(surface_coupling /= 0.0_real64) ) &
        error stop 'internal-coordinate driver coupling initialization failed'

    if( status_success /= 0 .or. status_derivative_unavailable /= 1 &
        .or. status_invalid_input /= 2 ) error stop 'status contract changed'
    if( primary_term_count /= 418 .or. secondary_term_count /= 330 ) &
        error stop 'H3 basis term counts changed'
    call get_h3_basis_term(1, first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [0,1,1,1]) ) &
        error stop 'first primary H3 basis term changed'
    call get_h3_basis_term(primary_term_count, first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [4,4,4,89]) ) &
        error stop 'last primary H3 basis term changed'
    call get_h3_basis_term(primary_term_count+1, first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [0,1,1,1]) ) &
        error stop 'first secondary H3 basis term changed'
    call get_h3_basis_term(primary_term_count+secondary_term_count, &
        first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [4,4,3,71]) ) &
        error stop 'last secondary H3 basis term changed'
    contains
    subroutine test_triatomic_kernel( distances_bohr, derivative_order, energy_hartree, gradient_hartree_per_bohr )
        real(real64), intent(in)  :: distances_bohr(3)
        integer, intent(in)       :: derivative_order
        real(real64), intent(out) :: energy_hartree, gradient_hartree_per_bohr(3)
        energy_hartree = sum(distances_bohr)
        gradient_hartree_per_bohr = 0.0_real64
        if( derivative_order == 1 ) gradient_hartree_per_bohr = [1.0_real64, 2.0_real64, 3.0_real64]
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
