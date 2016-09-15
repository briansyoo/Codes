!CASANDRA POST ANALYSIS LIBRARY
!FILE: cas_palib.f90


SUBROUTINE CAS_DENSITY(xyzfile, nslices, nspecies, nframes, begin_frame, end_frame, &
						Lx, Ly, Lz, &
						atype_mass, max_natoms, &
						natoms, nmolecules_array, density)
!**************************************************************************
!Routine to obtain the density profile along a cartesian axis
! Input variables: 
!		INTEGERS:
!		xyzfile - coordinate file containing all frames
!		nslices - number of slices to bin along the cartesian axis
!		nspecies - number of species contained in the system
!		nframes - number of frames we are analyzing
!		begin_frame - beginning frame index
!		end_frame - ending_frame index
!		REALS:
!		binwidth
!		ARRAYS:
!		max_natoms
!		L - array containing the length of the box for each frame
!		nmolecules_array - array containing number of molecules for each 
!		species in each frame
!		density (output) - density profile array
!
!**************************************************************************
IMPLICIT NONE
CHARACTER(120), INTENT(IN) :: xyzfile
INTEGER, INTENT(IN) :: nslices, nspecies, nframes, max_natoms, begin_frame, end_frame
REAL, INTENT(IN), DIMENSION(0:nframes-1) :: Lx, Ly, Lz
!f2py depend(nframes) :: Lx, Ly, Lz
REAL,INTENT(IN), DIMENSION(0:max_natoms-1,0:nspecies-1) :: atype_mass
!f2py depend(max_natoms,nspecies) :: atype_mass
INTEGER, INTENT(IN), DIMENSION(0:nspecies-1) :: natoms
!f2py depend(nspecies) :: natoms
INTEGER, INTENT(IN), DIMENSION(0:nframes-1,0:nspecies-1) :: nmolecules_array
!f2py depend(nframes,nspecies) :: nmolecules_array
REAL, INTENT(INOUT), DIMENSION(0:nslices-1) :: density
!f2py intent(in,out) :: density
!f2py depend(nslices) :: density
INTEGER :: i,j,k, nmolecules, slice
REAL :: rxp,ryp,rzp, binwidth, binvol,conv
CHARACTER(6) :: atype

conv = 1660.54!amu/Angstrom^3 to kg/m^3

OPEN(unit=5,file=xyzfile)

!loop over #frames
DO i=begin_frame,end_frame-1
	!first 2 lines is molecule #, and newline
	READ(*,*)
	READ(*,*)
	!define binwidth in the frame
	binwidth = REAL(Lz(i)/nslices)
	binvol = REAL(nslices/(Lx(i)*Ly(i)*Lz(i)))
	!loop over number of species
	DO j=0,nspecies-1
		nmolecules= nmolecules_array(i,j)*natoms(j)
		!loop over number of molecules
		DO k=0,nmolecules-1
			READ (*,*) atype, rxp,ryp,rzp
			slice = ABS(NINT((rzp+Lz(i)/2.0)/binwidth))
			density(slice) = density(slice)+1.0*binvol *atype_mass(MOD(k,natoms(j)),j)
		END DO

	END DO

END DO
CLOSE(unit = 5)
density = density/(end_frame-begin_frame)*conv

END SUBROUTINE





SUBROUTINE CAS_RDF
!Routine to generate the radial distribution function for a given index
IMPLICIT NONE
END SUBROUTINE


