!> @brief Evaluate the stateless EAM/DIM and modified-J1 CH4/Ni(111) energy kernel.
module chempotpy_ch4ni_runtime
    use chempotpy_kinds, only: wp
    implicit none
    private
    public :: evaluate_ch4ni_energy
    integer, parameter  :: nickel_type              = 1, hydrogen_type = 2, carbon_type = 3
    real(wp), parameter :: inner_cutoff_angstrom    = 5.0_wp, outer_cutoff_angstrom = 10.0_wp
    real(wp), parameter :: electrostatic_conversion = 14.3888_wp
    ! Preserve the released separated CH4+Ni(111) zero after adopting the shared CODATA conversions and PI.
    real(wp), parameter :: fixed_slab_zero_correction_hartree = 7.050451653896061e-4_wp
    integer, parameter  :: nickel_principal_numbers(10)       = [1, 1, 2, 2, 3, 3, 4, 4, 3, 3]
    integer, parameter  :: carbon_principal_numbers(4)        = [1, 1, 2, 2]
    real(wp), parameter :: nickel_factorials(10)              = [ &
        & 2.0_wp, 2.0_wp, 24.0_wp, 24.0_wp, 720.0_wp, 720.0_wp, 40320.0_wp, 40320.0_wp, 720.0_wp, 720.0_wp ]
    real(wp), parameter :: carbon_factorials(4)    = [2.0_wp, 2.0_wp, 24.0_wp, 24.0_wp]
    real(wp), parameter :: nickel_coefficients(10) = [ &
        & -0.00389_wp, -0.02991_wp, -0.03189_wp, 0.15289_wp, -0.20048_wp, -0.05423_wp, 0.49292_wp, &
        & 0.61875_wp, 0.4212_wp, 0.70658_wp ]
    real(wp), parameter :: carbon_coefficients(4) = [-0.20814_wp, -0.01071_wp, 0.08099_wp, 0.75045_wp]
    real(wp), parameter :: nickel_exponents(10)   = [ &
        & 54.88885_wp, 38.48431_wp, 27.42703_wp, 20.88204_wp, 10.95707_wp, 7.31958_wp, 3.9265_wp, &
        & 2.15289_wp, 12.67582_wp, 5.43253_wp ]
    real(wp), parameter :: carbon_exponents(4)            = [10.27273_wp, 17.91979_wp, 1.99841_wp, 2.88051_wp]
    real(wp), parameter :: nickel_radial_coefficients(10) = nickel_coefficients* &
        & (2.0_wp*nickel_exponents)**(real(nickel_principal_numbers, wp)+0.5_wp)/sqrt(nickel_factorials)
    real(wp), parameter :: carbon_radial_coefficients(4) = carbon_coefficients* &
        & (2.0_wp*carbon_exponents)**(real(carbon_principal_numbers, wp)+0.5_wp)/sqrt(carbon_factorials)
    real(wp), parameter :: hydrogen_radial_coefficient     = (2.0_wp*1.88972652_wp)**1.5_wp/sqrt(2.0_wp)
    real(wp), parameter :: nickel_embedding_parameters(14) = [ &
        & -303.28911_wp, 7.64744_wp, 87.98723_wp, 75.07472_wp, -522.77364_wp, 373.37826_wp, 0.1_wp, &
        & 0.11_wp, -2.014563533e10_wp, -538419456.0_wp, -4060125.0_wp, -11.03098054_wp, 39.42118_wp, &
        & 56.34512_wp ]

    contains

    !> @brief Return the CH4/Ni(111) energy in hartree for H,H,H,H,C Cartesian coordinates in angstrom.
    pure subroutine evaluate_ch4ni_energy( coordinates_angstrom, energy_hartree, components_hartree )
        use chempotpy_ch4ni_ni111_parameters, only: fixed_ni_coordinates_bohr
        use chempotpy_units,                  only: bohr_to_angstrom, hartree_to_electronvolt
        real(wp), intent(in)            :: coordinates_angstrom(5,3)
        real(wp), intent(out)           :: energy_hartree
        real(wp), intent(out), optional :: components_hartree(2)
        real(wp)                        :: bending_energy, dim_energy, nickel_embedding_energy
        real(wp)                        :: fixed_coordinates(size(fixed_ni_coordinates_bohr, 1),3)
        real(wp)                        :: fixed_density(size(fixed_ni_coordinates_bohr, 1)), gas_density(5)
        fixed_coordinates = fixed_ni_coordinates_bohr*bohr_to_angstrom
        call build_densities(coordinates_angstrom, fixed_coordinates, gas_density, fixed_density)
        call evaluate_eam_dim(coordinates_angstrom, fixed_coordinates, gas_density, fixed_density, dim_energy, &
            & nickel_embedding_energy)
        call evaluate_modified_j1(coordinates_angstrom/bohr_to_angstrom, gas_density, bending_energy)
        energy_hartree = (dim_energy+nickel_embedding_energy)/hartree_to_electronvolt+98.8029433_wp+ &
            & fixed_slab_zero_correction_hartree+bending_energy
        if( present(components_hartree) ) components_hartree = [ &
            & (dim_energy+nickel_embedding_energy)/hartree_to_electronvolt+98.8029433_wp+ &
            & fixed_slab_zero_correction_hartree, bending_energy ]
    end subroutine evaluate_ch4ni_energy

    !> @brief Accumulate the released directed atomic-density contributions at gas and fixed Ni sites.
    pure subroutine build_densities( gas_coordinates, fixed_coordinates, gas_density, fixed_density )
        real(wp), intent(in)  :: gas_coordinates(5,3), fixed_coordinates(:,:)
        real(wp), intent(out) :: gas_density(5), fixed_density(size(fixed_coordinates, 1))
        real(wp)              :: carbon_hydrogen_distance, distance, free_density_fraction
        integer               :: first_nickel, gas_atom, gas_type, second_nickel
        gas_density   = 0.0_wp
        fixed_density = 0.0_wp
        do first_nickel=1, size(fixed_coordinates, 1)-1
            do second_nickel=first_nickel+1, size(fixed_coordinates, 1)
                distance = norm2(fixed_coordinates(first_nickel,:)-fixed_coordinates(second_nickel,:))
                if( distance >= outer_cutoff_angstrom ) cycle
                fixed_density(first_nickel)  = fixed_density(first_nickel)+atomic_density(nickel_type, distance)
                fixed_density(second_nickel) = fixed_density(second_nickel)+atomic_density(nickel_type, distance)
            enddo
        enddo
        do gas_atom=1, 5
            gas_type = merge(hydrogen_type, carbon_type, gas_atom <= 4)
            if( gas_type == hydrogen_type ) &
                carbon_hydrogen_distance = norm2(gas_coordinates(gas_atom,:)-gas_coordinates(5,:))
            do first_nickel=1, size(fixed_coordinates, 1)
                distance = norm2(gas_coordinates(gas_atom,:)-fixed_coordinates(first_nickel,:))
                if( distance >= outer_cutoff_angstrom ) cycle
                gas_density(gas_atom) = gas_density(gas_atom)+atomic_density(nickel_type, distance)
                free_density_fraction = 1.0_wp
                if( gas_type == hydrogen_type ) free_density_fraction = available_hydrogen_fraction(distance, &
                    & carbon_hydrogen_distance)
                fixed_density(first_nickel) = fixed_density(first_nickel)+ free_density_fraction* &
                    & atomic_density(gas_type, distance)
            enddo
        enddo
    end subroutine build_densities

    !> @brief Evaluate the four generalized DIM roots and the nickel embedding sum in electronvolt.
    pure subroutine evaluate_eam_dim( gas_coordinates, fixed_coordinates, gas_density, fixed_density, &
            dim_energy, nickel_embedding_energy )
        use chempotpy_units, only: bohr_to_angstrom, hartree_to_electronvolt
        real(wp), intent(in)  :: gas_coordinates(5,3), fixed_coordinates(:,:), gas_density(5), fixed_density(:)
        real(wp), intent(out) :: dim_energy, nickel_embedding_energy
        real(wp)              :: exponential, gas_embedding(5), gas_singlet(5), gas_triplet(5)
        real(wp)              :: pair_sum, distance, sato_scale
        real(wp)              :: bond_singlet, bond_triplet, terms_singlet(3), terms_triplet(3)
        integer               :: gas_atom, gas_type, nickel
        do gas_atom=1, 4
            gas_embedding(gas_atom) = hydrogen_embedding(gas_density(gas_atom))
            if( gas_density(gas_atom) < 0.5_wp ) gas_embedding(gas_atom) = gas_embedding(gas_atom)+ &
                & symmetric_descriptor(-0.081682_wp, 1.0_wp, 5.0_wp, 0.0_wp)+ symmetric_descriptor(-0.006505_wp, &
                & 10.0_wp, 4.4_wp, 0.0_wp)
        enddo
        gas_embedding(5) = carbon_embedding(gas_density(5))
        do gas_atom=1, 5
            gas_type = merge(hydrogen_type, carbon_type, gas_atom <= 4)
            pair_sum = 0.0_wp
            do nickel=1, size(fixed_coordinates, 1)
                distance = norm2(gas_coordinates(gas_atom,:)-fixed_coordinates(nickel,:))
                if( distance < outer_cutoff_angstrom ) pair_sum = pair_sum+pair_potential(gas_type, nickel_type, distance)
            enddo
            gas_singlet(gas_atom) = gas_embedding(gas_atom)+pair_sum
            if( gas_type == hydrogen_type )then
                sato_scale = 0.5_wp*tanh(32.0_wp*(gas_density(gas_atom)-0.055_wp))+0.8_wp
            else
                sato_scale = 0.5_wp*tanh(64.2817_wp*(gas_density(gas_atom)-0.0069_wp))+0.5382_wp
            endif
            gas_triplet(gas_atom) = sato_scale*(-gas_embedding(gas_atom)+pair_sum)
        enddo
        dim_energy = 0.0_wp
        do gas_atom=1, 4
            distance     = norm2(gas_coordinates(gas_atom,:)-gas_coordinates(5,:))
            exponential  = exp(-(0.939_wp/bohr_to_angstrom)*(distance-2.0673_wp*bohr_to_angstrom))
            bond_singlet = 0.178852_wp*hartree_to_electronvolt*exponential*(exponential-2.0_wp)
            bond_triplet  = 0.5_wp*(1.0_wp-0.2699_wp)/(1.0_wp+0.2699_wp)*0.178852_wp* hartree_to_electronvolt* &
                & exponential*(exponential+2.0_wp)
            terms_singlet = [bond_singlet, gas_singlet(gas_atom), gas_singlet(5)]
            terms_triplet = [bond_triplet, gas_triplet(gas_atom), gas_triplet(5)]
            dim_energy    = dim_energy+lowest_dim_root(terms_singlet, terms_triplet)
        enddo
        nickel_embedding_energy = 0.0_wp
        do nickel=1, size(fixed_density)
            nickel_embedding_energy = nickel_embedding_energy+nickel_embedding(fixed_density(nickel))
        enddo
    end subroutine evaluate_eam_dim

    !> @brief Evaluate the modified J1 methane bending energy in hartree.
    pure subroutine evaluate_modified_j1( coordinates_bohr, gas_density, energy )
        real(wp), intent(in)  :: coordinates_bohr(5,3), gas_density(5)
        real(wp), intent(out) :: energy
        real(wp)              :: carbon_surface_distance, combination_energy, total_weight, weight
        real(wp)              :: hydrogen_surface_distances(4), weighted_energy, all_weights(20), all_energies(20)
        real(wp)              :: dummy_surface(3), entity_coordinates(6,3)
        integer               :: first_entity, combination, second_entity
        call density_pseudodistances(gas_density, hydrogen_surface_distances, carbon_surface_distance)
        entity_coordinates(1:5,:) = coordinates_bohr
        dummy_surface             = [coordinates_bohr(5,1), coordinates_bohr(5,2), 0.0_wp]
        entity_coordinates(6,:)   = dummy_surface
        combination               = 0
        do first_entity=1, 6
            if( first_entity == 5 ) cycle
            do second_entity=1, 6
                if( second_entity == 5 .or. second_entity == first_entity ) cycle
                combination = combination+1
                call evaluate_j1_combination(entity_coordinates, hydrogen_surface_distances, carbon_surface_distance, &
                    & first_entity, second_entity, combination_energy, weight)
                all_energies(combination) = combination_energy
                all_weights(combination)  = weight
            enddo
        enddo
        total_weight    = sum(all_weights)
        weighted_energy = 0.0_wp; weight = 0.0_wp
        do combination=1, 20
            if( all_weights(combination)/total_weight < 1.0e-8_wp ) cycle
            weighted_energy = weighted_energy+all_energies(combination)*all_weights(combination)
            weight          = weight+all_weights(combination)
        enddo
        energy = weighted_energy/weight
    end subroutine evaluate_modified_j1

    !> @brief Evaluate one ordered modified-J1 surface/atom assignment and its normalization weight.
    pure subroutine evaluate_j1_combination( coordinates, hydrogen_surface_distances, carbon_surface_distance, &
            first_entity, second_entity, energy, weight )
        use chempotpy_units, only: hartree_to_electronvolt
        real(wp), intent(in)  :: coordinates(6,3), hydrogen_surface_distances(4), carbon_surface_distance
        integer, intent(in)   :: first_entity, second_entity
        real(wp), intent(out) :: energy, weight
        integer               :: entities(6), candidate, count, index
        real(wp)              :: angles(6), angle_difference, angle_energy, attenuation(4), equilibrium_back
        real(wp)              :: equilibrium_front, force_constants(6), distances(15), squared_distances(15)
        real(wp)              :: exponential_first, exponential_second, force_scale, mean_distance
        entities(1) = first_entity; count = 1
        do candidate=1, 6
            if( candidate == 5 .or. candidate == first_entity .or. candidate == second_entity ) cycle
            count = count+1; entities(count) = candidate
        enddo
        entities(5) = 5; entities(6) = second_entity
        call j1_distances(coordinates, hydrogen_surface_distances, carbon_surface_distance, entities, distances)
        squared_distances = distances*distances
        call j1_angles(distances, squared_distances, angles)
        do index=1, 4
            exponential_first = 0.0_wp; exponential_second = 0.0_wp
            if( squared_distances(5+index) < 2500.0_wp ) &
                exponential_first = exp(-0.8999975_wp*squared_distances(5+index))
            if( distances(5+index)-1.402_wp < 50.0_wp ) &
                exponential_second = exp(-3.239921_wp*(distances(5+index)-1.402_wp)**2)
            attenuation(index) = (1.0_wp-exponential_first)*exp(-(0.448035_wp+ 0.606528_wp*exponential_second)* &
                & (distances(index)-2.0673_wp)**2)
        enddo
        force_scale = 2.55_wp+0.8_wp*exp(-3.0_wp*(distances(1)-2.0673_wp)**2)
        force_constants = force_scale*[attenuation(1)*attenuation(2), attenuation(1)*attenuation(3), attenuation(1)* &
            & attenuation(4), attenuation(2)*attenuation(3), attenuation(2)*attenuation(4), attenuation(3)* &
            & attenuation(4)]
        call equilibrium_angles(distances(1:4), equilibrium_front, equilibrium_back)
        angle_energy = 0.0_wp
        do index=1, 6
            angle_difference = angles(index)-merge(equilibrium_front, equilibrium_back, index <= 3)
            angle_energy     = angle_energy+0.5_wp*force_constants(index)*(angle_difference**2+angle_difference**6)
        enddo
        energy        = angle_energy/hartree_to_electronvolt+out_of_plane_energy(coordinates, entities, distances(1))
        mean_distance = sum(distances(10:12))/3.0_wp
        weight        = (mean_distance**6*(distances(5)**2/distances(6))**6)**2
    end subroutine evaluate_j1_combination

    !> @brief Construct the fifteen J1 pair distances for one ordered entity assignment.
    pure subroutine j1_distances( coordinates, hydrogen_surface_distances, carbon_surface_distance, entities, distances )
        real(wp), intent(in)  :: coordinates(6,3), hydrogen_surface_distances(4), carbon_surface_distance
        integer, intent(in)   :: entities(6)
        real(wp), intent(out) :: distances(15)
        integer, parameter    :: pair_indices(2,6) = reshape([1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])
        integer               :: index
        do index=1, 4
            distances(index) = entity_distance(coordinates, hydrogen_surface_distances, carbon_surface_distance, &
                & entities(index), 5)
        enddo
        distances(5) = entity_distance(coordinates, hydrogen_surface_distances, carbon_surface_distance, 5, entities(6))
        do index=1, 4
            distances(5+index) = entity_distance(coordinates, hydrogen_surface_distances, carbon_surface_distance, &
                & entities(index), entities(6))
        enddo
        do index=1, 6
            distances(9+index) = entity_distance(coordinates, hydrogen_surface_distances, carbon_surface_distance, &
                & entities(pair_indices(1,index)), entities(pair_indices(2,index)))
        enddo
    end subroutine j1_distances

    !> @brief Return an atom/atom or density-mapped atom/surface distance in bohr.
    pure real(wp) function entity_distance( coordinates, hydrogen_surface_distances, carbon_surface_distance, &
            first_entity, second_entity ) result(distance)
        real(wp), intent(in) :: coordinates(6,3), hydrogen_surface_distances(4), carbon_surface_distance
        integer, intent(in)  :: first_entity, second_entity
        integer              :: atom
        if( first_entity == 6 .or. second_entity == 6 )then
            atom = merge(second_entity, first_entity, first_entity == 6)
            if( atom == 5 )then
                distance = carbon_surface_distance
            else
                distance = hydrogen_surface_distances(atom)
            endif
        else
            distance = norm2(coordinates(first_entity,:)-coordinates(second_entity,:))
        endif
    end function entity_distance

    !> @brief Convert gas-site EAM densities into the five released surface pseudodistances.
    pure subroutine density_pseudodistances( density, hydrogen_distances, carbon_distance )
        real(wp), intent(in)  :: density(5)
        real(wp), intent(out) :: hydrogen_distances(4), carbon_distance
        real(wp)              :: carbon_a, carbon_b, carbon_c, carbon_weights(3)
        real(wp)              :: hydrogen_a, hydrogen_c, hydrogen_weights(3), reference_density
        integer               :: atom
        if( density(5) > 0.0_wp )then
            carbon_a = 1.4587_wp/(0.09393_wp**(-0.1551_wp)-0.09393_wp**0.0675_wp)* (density(5)**(-0.1551_wp)- &
                & density(5)**0.0675_wp)
            carbon_b = 2.067_wp/(0.023680_wp**(-0.4368_wp)-0.023680_wp**(-0.1816_wp))* (density(5)**(-0.4368_wp)- &
                & density(5)**(-0.1816_wp))
            carbon_c        = 1.4_wp
            carbon_weights  = [1.0_wp/density(5)**2, 4.0_wp/(density(5)-0.02385_wp)**2, 1.0_wp/(density(5)-0.18_wp)**2]
            carbon_distance = dot_product(carbon_weights, [carbon_a, carbon_b, carbon_c])/sum(carbon_weights)
        else
            carbon_distance = 200.0_wp
        endif
        do atom=1, 4
            if( density(atom) > 0.0_wp )then
                hydrogen_a = 1.4587_wp/(0.09393_wp**(-0.1551_wp)-0.09393_wp**0.0675_wp)* (density(atom)**(-0.1551_wp)- &
                    & density(atom)**0.0675_wp)
                hydrogen_c        = 1.4_wp
                reference_density = merge(0.010_wp, 0.01051_wp, atom == 1)
                hydrogen_weights = [1.0_wp/density(atom)**2, 8.0_wp/(density(atom)-reference_density)**2, 1.0_wp/ &
                    & (density(atom)-0.18_wp)**2]
                hydrogen_distances(atom) = dot_product(hydrogen_weights, [hydrogen_a, hydrogen_a, hydrogen_c])/ &
                    & sum(hydrogen_weights)
            else
                hydrogen_distances(atom) = 200.0_wp
            endif
        enddo
    end subroutine density_pseudodistances

    !> @brief Evaluate the six H-C-H angles from the J1 distance set.
    pure subroutine j1_angles( distances, squared_distances, angles )
        real(wp), intent(in)  :: distances(15), squared_distances(15)
        real(wp), intent(out) :: angles(6)
        real(wp)              :: cosine
        integer               :: index
        do index=1, 3
            cosine = (squared_distances(1)+squared_distances(index+1)-squared_distances(index+9))/ (2.0_wp* &
                & distances(1)*distances(index+1))
            angles(index) = acos(max(-0.999999999_wp, min(0.999999999_wp, cosine)))
        enddo
        do index=1, 2
            cosine = (squared_distances(2)+squared_distances(index+2)-squared_distances(index+12))/ (2.0_wp* &
                & distances(2)*distances(index+2))
            angles(index+3) = acos(max(-0.999999999_wp, min(0.999999999_wp, cosine)))
        enddo
        cosine    = (squared_distances(3)+squared_distances(4)-squared_distances(15))/ (2.0_wp*distances(3)* &
            & distances(4))
        angles(6) = acos(max(-0.999999999_wp, min(0.999999999_wp, cosine)))
    end subroutine j1_angles

    !> @brief Evaluate the modified J1 frontside and backside equilibrium angles.
    pure subroutine equilibrium_angles( carbon_hydrogen_distances, frontside, backside )
        use chempotpy_constants, only: PI
        real(wp), intent(in)  :: carbon_hydrogen_distances(4)
        real(wp), intent(out) :: frontside, backside
        real(wp), parameter   :: equilibrium_distance = 2.0673_wp, transition_distance = 6.275_wp
        real(wp), parameter   :: transition_width     = 0.1_wp, tetrahedral_angle = 1.9105588_wp
        real(wp)              :: a, b, c, difference, first_derivative, mean_distance, second_derivative
        mean_distance = sum(carbon_hydrogen_distances)/4.0_wp
        if( mean_distance > transition_distance )then
            frontside = PI/2.0_wp
        elseif( mean_distance <= equilibrium_distance )then
            frontside = tetrahedral_angle
        elseif( mean_distance <= equilibrium_distance+transition_width )then
            first_derivative = -4.0_wp*0.080768_wp/transition_width**2
            a                = 3.0_wp*first_derivative/transition_width**2
            b                = -(5.0_wp*a*transition_width**2+first_derivative)/(2.0_wp*transition_width)
            c                = -(20.0_wp*a*transition_width**2+12.0_wp*b*transition_width)/6.0_wp
            difference       = mean_distance-equilibrium_distance
            frontside        = a*difference**5+b*difference**4+c*difference**3+tetrahedral_angle
        elseif( mean_distance <= transition_distance-transition_width )then
            c         = 1.0133_wp*(mean_distance-equilibrium_distance)*exp(-7.9866e-4_wp*mean_distance**3)
            frontside = tetrahedral_angle+(tetrahedral_angle-PI/2.0_wp)*((1.0_wp-tanh(c))-1.0_wp)
        else
            first_derivative = (tetrahedral_angle-4.0_wp*0.080768_wp*(transition_distance-transition_width- &
                & equilibrium_distance)-PI/2.0_wp)/transition_width**3
            second_derivative = -4.0_wp*0.080768_wp/transition_width**2
            a                 = -3.0_wp*(2.0_wp*first_derivative+second_derivative)/transition_width**2
            b                 = (second_derivative+5.0_wp*a*transition_width**2)/(2.0_wp*transition_width)
            c                 = (12.0_wp*b*transition_width-20.0_wp*a*transition_width**2)/6.0_wp
            difference        = mean_distance-transition_distance
            frontside         = a*difference**5+b*difference**4+c*difference**3+PI/2.0_wp
        endif
        a        = (cos(tetrahedral_angle)+0.5_wp)/cos(tetrahedral_angle)**2
        backside = acos(a*cos(frontside)**2-0.5_wp)
    end subroutine equilibrium_angles

    !> @brief Evaluate the modified J1 methyl out-of-plane term for one assignment.
    pure real(wp) function out_of_plane_energy( coordinates, entities, mapped_first_distance ) result(energy)
        use chempotpy_constants, only: PI
        real(wp), intent(in) :: coordinates(6,3), mapped_first_distance
        integer, intent(in)  :: entities(6)
        real(wp)             :: cross(3), delta, force_quadratic, force_quartic, gamma, normal
        real(wp)             :: plane_first(3), plane_second(3), phi_zero, radial_vector(3), rest
        real(wp)             :: switching, vector_norm
        integer              :: index
        plane_first  = coordinates(entities(3),:)-coordinates(entities(2),:)
        plane_second = coordinates(entities(4),:)-coordinates(entities(2),:)
        cross        = cross_product_3d(plane_first, plane_second)
        normal       = norm2(cross)
        rest         = -3.9568698e-5_wp*mapped_first_distance**3
        switching    = 0.0_wp
        if( rest < 70.0_wp ) switching = 1.0_wp-tanh(0.2509214_wp*(mapped_first_distance-2.052_wp)*exp(rest))
        phi_zero = 1.9106_wp+(1.9106_wp-PI/2.0_wp)*(switching-1.0_wp)
        switching      = 1.0_wp-tanh(0.02102956_wp*(mapped_first_distance-2.052_wp)* &
            & (mapped_first_distance-0.5798533_wp)**2)
        force_quadratic = (1.0_wp-switching)*0.0219_wp
        force_quartic   = (1.0_wp-switching)*0.0438_wp
        energy          = 0.0_wp
        do index=2, 4
            radial_vector = coordinates(entities(index),:)-coordinates(entities(5),:)
            vector_norm   = norm2(radial_vector)
            gamma         = dot_product(cross, radial_vector)/(normal*vector_norm)
            delta         = acos(max(-1.0_wp, min(1.0_wp, gamma)))-phi_zero
            energy        = energy+force_quadratic*delta**2+force_quartic*delta**4
        enddo
    end function out_of_plane_energy

    !> @brief Return the three-dimensional cross product used by the modified J1 kernel.
    pure function cross_product_3d( first, second ) result(product)
        real(wp), intent(in) :: first(3), second(3)
        real(wp)             :: product(3)
        product = [first(2)*second(3)-first(3)*second(2), first(3)*second(1)-first(1)*second(3), first(1)*second(2)- &
            & first(2)*second(1)]
    end function cross_product_3d

    !> @brief Return one released atomic-density contribution after the smooth radial cutoff.
    pure real(wp) function atomic_density( atom_type, distance ) result(density)
        use chempotpy_constants, only: PI
        integer, intent(in)  :: atom_type
        real(wp), intent(in) :: distance
        real(wp)             :: electron_count_d, electron_count_s, radial_d, radial_s
        integer              :: orbital
        radial_s = 0.0_wp; radial_d = 0.0_wp
        select case( atom_type )
        case( nickel_type )
            do orbital=1, 8
                radial_s = radial_s+nickel_radial_coefficients(orbital)*distance** (nickel_principal_numbers(orbital)- &
                    & 1)*exp(-nickel_exponents(orbital)*distance)
            enddo
            do orbital=9, 10
                radial_d = radial_d+nickel_radial_coefficients(orbital)*distance** (nickel_principal_numbers(orbital)- &
                    & 1)*exp(-nickel_exponents(orbital)*distance)
            enddo
            electron_count_s = 2.0_wp; electron_count_d = 8.0_wp
        case( hydrogen_type )
            radial_s         = hydrogen_radial_coefficient*exp(-1.88972652_wp*distance)
            electron_count_s = 1.0_wp; electron_count_d = 0.0_wp
        case( carbon_type )
            do orbital=1, 4
                radial_s = radial_s+carbon_radial_coefficients(orbital)*distance** (carbon_principal_numbers(orbital)- &
                    & 1)*exp(-carbon_exponents(orbital)*distance)
            enddo
            electron_count_s = 2.0_wp; electron_count_d = 2.0_wp
        end select
        density = (electron_count_s*radial_s**2+electron_count_d*radial_d**2)/(4.0_wp*PI)*cutoff(distance)
    end function atomic_density

    !> @brief Return the fraction of one H density available to a neighboring Ni atom.
    pure real(wp) function available_hydrogen_fraction( hydrogen_nickel_distance, carbon_hydrogen_distance ) result(fraction)
        real(wp), intent(in) :: hydrogen_nickel_distance, carbon_hydrogen_distance
        real(wp)             :: carbon_factor, nickel_factor
        nickel_factor = exp(-1.464_wp*(hydrogen_nickel_distance-1.9_wp))
        carbon_factor = exp(-1.774_wp*(carbon_hydrogen_distance-1.09_wp))
        fraction      = nickel_factor/(nickel_factor+carbon_factor)
    end function available_hydrogen_fraction

    !> @brief Return the released EAM7-form hydrogen embedding energy in electronvolt.
    pure real(wp) function hydrogen_embedding( density ) result(energy)
        real(wp), intent(in) :: density
        real(wp), parameter  :: boundary = 0.009_wp, width = 0.007_wp
        if( density <= boundary-width )then
            energy = merge(-4.8029_wp*density**0.38553_wp*exp(-34.462_wp*density), 0.0_wp, density > 0.0_wp)
        elseif( density <= boundary )then
            energy = 0.1840945234e12_wp*(density-boundary)**5+0.3276017408e10_wp*(density-boundary)**4+ 15102027.0_wp* &
                & (density-boundary)**3-0.1303885824e12_wp*(density-(boundary-width))**5+ 0.2375889152e10_wp*(density- &
                & (boundary-width))**4-11946214.0_wp*(density-(boundary-width))**3
        else
            energy = -74.556_wp*density**1.0154_wp*exp(-7.2784_wp*density)
        endif
    end function hydrogen_embedding

    !> @brief Return the released piecewise carbon embedding energy in electronvolt.
    pure real(wp) function carbon_embedding( density ) result(energy)
        real(wp), intent(in) :: density
        if( density <= 0.002_wp )then
            energy = merge(-6.4703_wp*density**0.41553_wp*exp(16.04_wp*density), 0.0_wp, density > 0.0_wp)
        elseif( density <= 0.003_wp )then
            energy = 0.3355955375e16_wp*(density-0.003_wp)**5+0.8340164313e13_wp*(density-0.003_wp)**4+ &
                & 0.5489281024e10_wp*(density-0.003_wp)**3-0.3223638975e16_wp*(density-0.002_wp)**5+ &
                & 0.8124186493e13_wp*(density-0.002_wp)**4-0.5501493760e10_wp*(density-0.002_wp)**3
        elseif( density <= 0.015_wp )then
            energy = -11.236_wp*density**0.52239_wp*exp(35.422_wp*density)
        elseif( density <= 0.022_wp )then
            energy = 0.9558951199e12_wp*(density-0.022_wp)**5+0.1648154112e11_wp*(density-0.022_wp)**4+ 74744608.0_wp* &
                & (density-0.022_wp)**3-0.9270320824e12_wp*(density-0.015_wp)**5+ 0.1628130099e11_wp*(density- &
                & 0.015_wp)**4-76418312.0_wp*(density-0.015_wp)**3
        else
            energy = -119.2417_wp*density**0.8279_wp*exp(-28.5356_wp*density)
        endif
    end function carbon_embedding

    !> @brief Return the EAM6/EAM7 nickel embedding energy in electronvolt.
    pure real(wp) function nickel_embedding( density ) result(energy)
        real(wp), intent(in) :: density
        real(wp)             :: delta, scaled_density
        if( density > nickel_embedding_parameters(8) )then
            energy = nickel_embedding_parameters(12)
        elseif( density > nickel_embedding_parameters(7) )then
            delta = density-nickel_embedding_parameters(8)
            energy = nickel_embedding_parameters(12)+delta**3*(nickel_embedding_parameters(11)+ delta* &
                & (nickel_embedding_parameters(10)+delta*nickel_embedding_parameters(9)))
        else
            scaled_density = nickel_embedding_parameters(13)*density
            energy = nickel_embedding_parameters(1)*exp(-nickel_embedding_parameters(2)*density)*density+ &
                & (nickel_embedding_parameters(3)*density)**3*exp(-nickel_embedding_parameters(4)*density)+ &
                & nickel_embedding_parameters(5)*exp(-nickel_embedding_parameters(6)*density)*density+ &
                & scaled_density**5*exp(-nickel_embedding_parameters(14)*density)
        endif
    end function nickel_embedding

    !> @brief Return one symmetric reciprocal-exponential descriptor energy.
    pure real(wp) function symmetric_descriptor( amplitude, slope, center, value ) result(energy)
        real(wp), intent(in) :: amplitude, slope, center, value
        energy = amplitude/(exp(slope*(value-center))+exp(slope*(center-value)))
    end function symmetric_descriptor

    !> @brief Return the lowest generalized eigenvalue of the released two-state DIM problem.
    pure real(wp) function lowest_dim_root( singlet, triplet ) result(root)
        real(wp), intent(in) :: singlet(3), triplet(3)
        real(wp)             :: coefficient_a, coefficient_b, coefficient_c, discriminant
        real(wp)             :: exchange(3), hamiltonian(2,2), mean(3), mean_sum, surface_exchange
        mean             = 0.5_wp*(singlet+triplet)
        exchange         = 0.5_wp*(singlet-triplet)
        mean_sum         = sum(mean)
        surface_exchange = exchange(2)+exchange(3)
        hamiltonian(1,1) = 2.0_wp*(2.0_wp*(mean_sum+exchange(1))-surface_exchange)
        hamiltonian(1,2) = -2.0_wp*(mean_sum+exchange(1)+surface_exchange)
        hamiltonian(2,1) = hamiltonian(1,2)
        hamiltonian(2,2) = 2.0_wp*(2.0_wp*(mean_sum+surface_exchange)-exchange(1))
        coefficient_a    = 12.0_wp
        coefficient_b    = -4.0_wp*(hamiltonian(1,1)+hamiltonian(1,2)+hamiltonian(2,2))
        coefficient_c    = hamiltonian(1,1)*hamiltonian(2,2)-hamiltonian(1,2)**2
        discriminant     = max(0.0_wp, coefficient_b**2-4.0_wp*coefficient_a*coefficient_c)
        root             = (-coefficient_b-sqrt(discriminant))/(2.0_wp*coefficient_a)
    end function lowest_dim_root

    !> @brief Return one screened pair interaction in electronvolt.
    pure real(wp) function pair_potential( first_type, second_type, distance ) result(energy)
        integer, intent(in)  :: first_type, second_type
        real(wp), intent(in) :: distance
        energy = electrostatic_conversion*effective_charge(first_type, distance)* effective_charge(second_type, &
            & distance)/distance
    end function pair_potential

    !> @brief Return the released Ni, H, or C effective charge.
    pure real(wp) function effective_charge( atom_type, distance ) result(charge)
        integer, intent(in)  :: atom_type
        real(wp), intent(in) :: distance
        select case( atom_type )
        case( nickel_type )
            charge = 10.0_wp*exp(-1.8633_wp*distance)*(1.0_wp+0.8957_wp*distance)*cutoff(distance)
        case( hydrogen_type )
            if( distance <= 1.09_wp )then
                charge = 1.0_wp-5.957_wp*distance**5+9.804_wp*distance**4-3.705_wp*distance**3
            else
                charge = (distance/0.93232_wp)**2.77511_wp*exp(-distance/0.37473_wp)+ (distance/0.54879_wp)** &
                    & 142.67761_wp*exp(-distance/0.01012_wp)+ (distance/0.37697_wp)**95.4552_wp*exp(-distance/ &
                    & 0.01073_wp)+ 0.015_wp/(exp(8.0_wp*(distance-2.4_wp))+exp(8.0_wp*(2.4_wp-distance)))
            endif
        case( carbon_type )
            charge = 26.2416_wp*(0.3736_wp*exp(-3.3380_wp*distance)- 5.0788_wp*distance*exp(-14.0996_wp*distance))
        end select
    end function effective_charge

    !> @brief Return the released smooth radial cutoff.
    pure real(wp) function cutoff( distance ) result(value)
        use chempotpy_constants, only: PI
        real(wp), intent(in) :: distance
        if( distance <= inner_cutoff_angstrom )then
            value = 1.0_wp
        elseif( distance < outer_cutoff_angstrom )then
            value = (distance-outer_cutoff_angstrom)/(inner_cutoff_angstrom-outer_cutoff_angstrom)+ sin(PI*((distance- &
                & inner_cutoff_angstrom)+(distance-outer_cutoff_angstrom))/ (inner_cutoff_angstrom- &
                & outer_cutoff_angstrom))/(2.0_wp*PI)
        else
            value = 0.0_wp
        endif
    end function cutoff
end module chempotpy_ch4ni_runtime
