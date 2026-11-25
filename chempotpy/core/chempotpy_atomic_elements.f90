!> @brief Provide canonical periodic-table symbols, atomic numbers, and standard atomic masses.
module chempotpy_atomic_elements
    use chempotpy_kinds, only: wp
    implicit none
    private
    integer, parameter, public  :: element_count                  = 118
    character(len=2), parameter :: element_symbols(element_count) = [ &
        & 'H ', 'He', &
        & 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
        & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
        & 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
        & 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
        & 'Cs', 'Ba', 'La', &
        &           'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
        &                'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
        & 'Fr', 'Ra', 'Ac', &
        &           'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', &
        &                'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' ]
    ! NIST standard atomic weights or representative isotope masses, in unified atomic mass units.
    real(wp), parameter :: standard_atomic_masses_u(element_count) = [ &
        &   1.00794075_wp, 4.00260193_wp, 6.94003660_wp, 9.01218307_wp, &
        &  10.81102805_wp, 12.01073590_wp, 14.00670321_wp, 15.99940492_wp, &
        &  18.99840316_wp, 20.18004638_wp, 22.98976928_wp, 24.30505162_wp, &
        &  26.98153853_wp, 28.08549871_wp, 30.97376200_wp, 32.06478741_wp, &
        &  35.45293758_wp, 39.94779856_wp, 39.09830091_wp, 40.07802251_wp, &
        &  44.95590828_wp, 47.86674496_wp, 50.94146504_wp, 51.99613176_wp, &
        &  54.93804391_wp, 55.84514443_wp, 58.93319429_wp, 58.69334711_wp, &
        &  63.54603995_wp, 65.37778253_wp, 69.72306607_wp, 72.62755016_wp, &
        &  74.92159457_wp, 78.95938856_wp, 79.90352778_wp, 83.79800000_wp, &
        &  85.46766360_wp, 87.61664447_wp, 88.90584030_wp, 91.22364160_wp, &
        &  92.90637300_wp, 95.95978854_wp, 97.90721240_wp, 101.06494014_wp, &
        & 102.90549800_wp, 106.41532751_wp, 107.86814963_wp, 112.41155782_wp, &
        & 114.81808663_wp, 118.71011259_wp, 121.75978367_wp, 127.60312648_wp, &
        & 126.90447190_wp, 131.29276145_wp, 132.90545196_wp, 137.32689163_wp, &
        & 138.90546887_wp, 140.11573074_wp, 140.90765760_wp, 144.24159603_wp, &
        & 144.91275590_wp, 150.36635571_wp, 151.96437813_wp, 157.25213065_wp, &
        & 158.92535470_wp, 162.49947282_wp, 164.93032880_wp, 167.25908265_wp, &
        & 168.93421790_wp, 173.05415017_wp, 174.96681496_wp, 178.48497872_wp, &
        & 180.94787564_wp, 183.84177755_wp, 186.20670455_wp, 190.22485963_wp, &
        & 192.21605165_wp, 195.08445686_wp, 196.96656879_wp, 200.59916703_wp, &
        & 204.38341284_wp, 207.21690806_wp, 208.98039910_wp, 208.98243080_wp, &
        & 209.98714790_wp, 222.01757820_wp, 223.01973600_wp, 226.02541030_wp, &
        & 227.02775230_wp, 232.03805580_wp, 231.03588420_wp, 238.02891046_wp, &
        & 237.04817360_wp, 244.06420530_wp, 243.06138130_wp, 247.07035410_wp, &
        & 247.07030730_wp, 251.07958860_wp, 252.08298000_wp, 257.09510610_wp, &
        & 258.09843150_wp, 259.10103000_wp, 262.10961000_wp, 267.12179000_wp, &
        & 269.12791000_wp, 271.13393000_wp, 270.13336000_wp, 276.14846000_wp, &
        & 276.15159000_wp, 280.16131000_wp, 282.16912000_wp, 284.17416000_wp, &
        & 284.17873000_wp, 289.19042000_wp, 288.19274000_wp, 293.20449000_wp, &
        & 292.20746000_wp, 294.21392000_wp ]
    public :: all_atomic_numbers_valid
    public :: atomic_mass_u
    public :: atomic_number
    public :: atomic_symbol
    public :: is_valid_atomic_number
    public :: is_valid_atomic_symbol
    interface atomic_mass_u
        module procedure atomic_mass_u_from_number
        module procedure atomic_mass_u_from_symbol
    end interface atomic_mass_u

    contains

    !> @brief Return the atomic number for a case-insensitive element symbol, or zero when unknown.
    pure function atomic_number( symbol ) result(number)
        character(len=*), intent(in) :: symbol
        character(len=2)             :: normalized_symbol
        integer                      :: number
        integer                      :: index
        normalized_symbol = lowercase_symbol(symbol)
        number            = 0
        do index=1, element_count
            if( normalized_symbol == lowercase_symbol(element_symbols(index)) )then
                number = index
                return
            endif
        enddo
        if( normalized_symbol == 'd ' .or. normalized_symbol == 't ' ) number = 1
    end function atomic_number

    !> @brief Return the canonical symbol for an atomic number, or two blanks when invalid.
    pure elemental function atomic_symbol( number ) result(symbol)
        integer, intent(in) :: number
        character(len=2)    :: symbol
        symbol = '  '
        if( is_valid_atomic_number(number) ) symbol = element_symbols(number)
    end function atomic_symbol

    !> @brief Return whether an atomic number identifies one of the 118 recognized elements.
    pure elemental logical function is_valid_atomic_number( number )
        integer, intent(in) :: number
        is_valid_atomic_number = number >= 1 .and. number <= element_count
    end function is_valid_atomic_number

    !> @brief Return whether a symbol identifies an element; D and T are accepted hydrogen aliases.
    pure logical function is_valid_atomic_symbol( symbol )
        character(len=*), intent(in) :: symbol
        is_valid_atomic_symbol = atomic_number(symbol) /= 0
    end function is_valid_atomic_symbol

    !> @brief Return whether every value in an array is a recognized atomic number.
    pure logical function all_atomic_numbers_valid( numbers )
        integer, intent(in) :: numbers(:)
        all_atomic_numbers_valid = all(is_valid_atomic_number(numbers))
    end function all_atomic_numbers_valid

    !> @brief Return the standard atomic mass in unified atomic mass units, or zero when invalid.
    pure elemental function atomic_mass_u_from_number( number ) result(mass)
        integer, intent(in) :: number
        real(wp)            :: mass
        mass = 0.0_wp
        if( is_valid_atomic_number(number) ) mass = standard_atomic_masses_u(number)
    end function atomic_mass_u_from_number

    !> @brief Return the standard atomic mass in unified atomic mass units, or zero when unknown.
    pure function atomic_mass_u_from_symbol( symbol ) result(mass)
        character(len=*), intent(in) :: symbol
        real(wp)                     :: mass
        mass = atomic_mass_u_from_number(atomic_number(symbol))
    end function atomic_mass_u_from_symbol

    !> @brief Normalize a one- or two-letter element symbol to lowercase ASCII.
    pure function lowercase_symbol( input_symbol ) result(output_symbol)
        character(len=*), intent(in)  :: input_symbol
        character(len=2)              :: output_symbol
        character(len=:), allocatable :: adjusted_symbol
        integer                       :: character_code
        integer                       :: index
        output_symbol   = '  '
        adjusted_symbol = trim(adjustl(input_symbol))
        if( len(adjusted_symbol) < 1 .or. len(adjusted_symbol) > 2 ) return
        do index=1, len(adjusted_symbol)
            character_code = iachar(adjusted_symbol(index:index))
            if( character_code >= iachar('A') .and. character_code <= iachar('Z') ) &
                character_code = character_code+iachar('a')-iachar('A')
            if( character_code < iachar('a') .or. character_code > iachar('z') ) return
            output_symbol(index:index) = achar(character_code)
        enddo
    end function lowercase_symbol
end module chempotpy_atomic_elements
