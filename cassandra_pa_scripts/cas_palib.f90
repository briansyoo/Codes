!CASANDRA POST ANALYSIS LIBRARY
!FILE: cas_palib.f90
!
!**************************************************************************
! This library contains all of the analysis codes called by each python
! script. Make sure to compile the code using the install script.
! This library is implemented using the python wrapper F2PY. 
! 
! Writte by: Brian Yoo
!**************************************************************************

SUBROUTINE CAS_DENSITY(xyzfile, nslices, nspecies, &
						nframes, begin_frame, end_frame, &
						binwidth, &
						Lx, Ly, Lz, dime, &
						atype_mass, max_natoms, &
						natoms, nmolecules, density)
!**************************************************************************
!Routine to obtain the density profile along a cartesian axis
!       Input Variables: 
!
!		xyzfile (INT) - coordinate file containing all frames
!		nslices (INT) - number of slices to bin along the cartesian axis
!		nspecies (INT) - number of species contained in the system
!		nframes (INT) - number of frames we are analyzing
!		begin_frame (INT) - beginning frame index
!		end_frame (INT) - ending_frame index
!		Lx, Ly, Lz (ARRAY) - array containing the length of the box for each frame
!		dime (CHAR) - dimension of axis we will perform analysis on
!		atype_mass (MATRIX) - matrix containing masses of atoms for each species
!		max_natoms (INT) - max number of atoms out of all species
!		natoms (ARRAY) - array containing number of atoms for each species
!		nmolecules (MATRIX) - array containing number of molecules for each 
!		species in each frame
!		density (ARRAY) - density profile array
!
!**************************************************************************
IMPLICIT NONE
CHARACTER(120), INTENT(IN) :: xyzfile, dime
INTEGER, INTENT(IN) :: nslices, nspecies, nframes, max_natoms, begin_frame, end_frame
REAL, INTENT(IN) :: binwidth
REAL, INTENT(IN), DIMENSION(0:nframes-1) :: Lx, Ly, Lz
!f2py depend(nframes) :: Lx, Ly, Lz
REAL,INTENT(IN), DIMENSION(0:max_natoms-1,0:nspecies-1) :: atype_mass
!f2py depend(max_natoms,nspecies) :: atype_mass
INTEGER, INTENT(IN), DIMENSION(0:nspecies-1) :: natoms
!f2py depend(nspecies) :: natoms
INTEGER, INTENT(IN), DIMENSION(0:nframes-1,0:nspecies-1) :: nmolecules
!f2py depend(nframes,nspecies) :: nmolecules
REAL, INTENT(INOUT), DIMENSION(0:nslices) :: density !added extra slice
!f2py intent(in,out) :: density
!f2py depend(nslices) :: density

INTEGER :: i,j,k,l, slice
REAL :: rxp,ryp,rzp, binvol,conv
REAL, ALLOCATABLE, DIMENSION(:) :: Lout
CHARACTER(6) :: atype

conv = 1660.54!amu/Angstrom^3 to kg/m^3 conversion factor

OPEN(unit=5,file=xyzfile)

ALLOCATE(Lout(0:nframes-1))


!loop over #frames
DO i=0,end_frame-1
	!first 2 lines is molecule #, and newline
	READ(*,*)
	READ(*,*)

	!loop over number of species
	DO j=0,nspecies-1
		!loop over number of molecules
		DO k=0,nmolecules(i,j)-1
			DO l=0, natoms(j)-1 
				READ (*,*) atype, rxp,ryp,rzp
				IF (i>=begin_frame) THEN
					binvol = REAL(nslices/(Lx(i)*Ly(i)*Lz(i)))

					!wrap coordinates
					IF (rxp < -Lx(i)/2.0) THEN
						rxp = rxp + Lx(i)
					ELSEIF (rxp > Lx(i)/2.0) THEN
						rxp = rxp - Lx(i)
					ENDIF
					IF (ryp < -Ly(i)/2.0) THEN
						ryp = ryp + Ly(i)
					ELSEIF (ryp > Ly(i)/2.0) THEN
						ryp = ryp - Ly(i)
					END IF
					IF (rzp < -Lz(i)/2.0) THEN
						rzp = rzp + Lz(i)
					ELSEIF (rzp > Lz(i)/2.0) THEN
						rzp = rzp - Lz(i)
					END IF 
					
					!determine slice index					
					IF (dime == 'Z' .OR. dime == 'z') THEN
						slice = NINT((rzp+REAL(Lz(i)/2.0))/binwidth)
					ELSEIF (dime == 'Y' .OR. dime =='y') THEN
						slice = NINT((ryp+REAL(Ly(i)/2.0))/binwidth)
					ELSEIF (dime == 'X' .OR. dime =='x') THEN
						slice = NINT((rxp+REAL(Lx(i)/2.0))/binwidth)
					END IF


					density(slice) = density(slice)+1.0*binvol *atype_mass(l,j)
				END IF
			END DO
		END DO

	END DO

END DO

CLOSE(unit = 5)
density = density/(end_frame-begin_frame)*conv


END SUBROUTINE



SUBROUTINE CAS_RDF(xyzfile,nbins,binwidth,nspecies,&
					nframes,begin_frame,end_frame,&
					Lx,Ly,Lz,Lmin,&
					species1_ndx, species2_ndx, &
					atype1_ndx, atype2_ndx, &
					atype1_name, atype2_name, &
					atype_mass,max_natoms,max_nmols,natoms,&
					nmolecules,com_flag, rdf)
!**************************************************************************
!Routine to generate the radial distribution function for given pair indices
!
! 		Input Variables: 
!		xyzfile (INT) - coordinate file containing all frames
!		nbins (INT) - number of bins
!		binwidth (REAL) - based on minimum box length
!		nspecies (INT) - number of species contained in the system
!		nframes (INT) - number of frames we are analyzing
!		begin_frame (INT) - beginning frame index
!		end_frame (INT) - ending_frame index
!		Lx, Ly, Lz (ARRAY) - box lengths 
!		Lmin (REAL) - global minimum box length from Lx,Ly, and Lz
!		species1_ndx, species2_ndx (INT) - species index
!		atype1_name, atype_name (CHAR) -name of atomtype
!		atype_mass (MATRIX) - contains masses of each atomtype for each species
!		max_natoms (INT) - max number of atoms out of atoms in species
!		max_nmols (INT) - max number of molecules out of molecules in species
!		natoms (ARRAY) - number of atoms in each species
!		nmolecules (MATRIX) - array containing number of molecules for each 
!		species in each frame
!		com_flag (LOGICAL) - flag determining analysis: com-com RDF or atom-atom RDF
!		rdf (output) - rdf array
!
!**************************************************************************
IMPLICIT NONE
CHARACTER(120), INTENT(IN) :: xyzfile
INTEGER, INTENT(IN) :: nbins, nspecies, nframes, begin_frame, end_frame, &
						species1_ndx, species2_ndx, atype1_ndx, atype2_ndx, max_natoms, max_nmols
CHARACTER(6), INTENT(IN) :: atype1_name, atype2_name
REAL, INTENT(IN) :: binwidth, Lmin
REAL, INTENT(IN), DIMENSION(0:nframes-1) :: Lx, Ly, Lz
!f2py depend(nframes) :: Lx, Ly, Lz
REAL,INTENT(IN), DIMENSION(0:max_natoms-1,0:nspecies-1) :: atype_mass
!f2py depend(max_natoms,nspecies) :: atype_mass
INTEGER, INTENT(IN), DIMENSION(0:nspecies-1) :: natoms
!f2py depend(nspecies) :: natoms
INTEGER, INTENT(IN), DIMENSION(0:nframes-1,0:nspecies-1) :: nmolecules
!f2py depend(nframes,nspecies) :: nmolecules
LOGICAL, INTENT(IN) :: com_flag
REAL, INTENT(INOUT), DIMENSION(0:nbins-1) :: rdf
!f2py intent(in,out) :: rdf
!f2py depend(nbins) :: rdf


INTEGER :: i, j, k, l, m, n, p, ia, ja
INTEGER :: bin_ndx
REAL :: rxp, ryp, rzp, rsq, shell_volume, bulk_density
REAL :: rxij, ryij, rzij, rxijp, ryijp, rzijp
CHARACTER(6) :: atype
REAL, ALLOCATABLE, DIMENSION(:) :: atype1_rx, atype1_ry, atype1_rz, atype2_rx, atype2_ry, atype2_rz
REAL, ALLOCATABLE, DIMENSION(:) :: mtype1_rxcom, mtype1_rycom, mtype1_rzcom, &
									mtype2_rxcom, mtype2_rycom, mtype2_rzcom
REAL :: mtype1_MW, mtype2_MW
REAL, ALLOCATABLE, DIMENSION(:) :: ncounts



ALLOCATE(ncounts(0:nbins-1))


bulk_density = 0.0
ncounts(:) = 0.0
rdf(:) = 0.0
OPEN(unit=5,file=xyzfile)

!CALCULATE CENTER-OF-MASS RDF
IF (com_flag .eqv. .TRUE.) THEN
	ALLOCATE(mtype1_rxcom(0:max_nmols-1),mtype1_rycom(0:max_nmols-1),mtype1_rzcom(0:max_nmols-1), &
		mtype2_rxcom(0:max_nmols-1),mtype2_rycom(0:max_nmols-1),mtype2_rzcom(0:max_nmols-1))


	mtype1_MW = SUM(atype_mass(:,species1_ndx))
	mtype2_MW = SUM(atype_mass(:,species2_ndx))
	!read xyz file lines for each frame
	DO i=0,end_frame-1
		mtype1_rxcom(:) = 0.0
		mtype1_rycom(:) = 0.0
		mtype1_rzcom(:) = 0.0
		mtype2_rxcom(:) = 0.0
		mtype2_rycom(:) = 0.0
		mtype2_rzcom(:) = 0.0
		!reset counter for atype index
		!first 2 lines is molecule #, and newline
		READ(*,*)
		READ(*,*)
		!loop over number of species
		DO j=0,nspecies-1
			!loop over number of molecules to store coordinates of atomtype
			DO k=0,nmolecules(i,j)-1
				DO l=0, natoms(j)-1
					READ (*,*) atype, rxp,ryp,rzp
					!calculate center of mass for each molecule
					IF (species1_ndx == j) THEN
						mtype1_rxcom(k) = mtype1_rxcom(k) + atype_mass(l,j)*rxp
						mtype1_rycom(k) = mtype1_rycom(k) + atype_mass(l,j)*ryp
						mtype1_rzcom(k) = mtype1_rzcom(k) + atype_mass(l,j)*rzp
					ELSEIF (species2_ndx == j) THEN
						mtype2_rxcom(k) = mtype2_rxcom(k) + atype_mass(l,j)*rxp
						mtype2_rycom(k) = mtype2_rycom(k) + atype_mass(l,j)*ryp
						mtype2_rzcom(k) = mtype2_rzcom(k) + atype_mass(l,j)*rzp
					ENDIF
				END DO

			END DO
		
		END DO

		!if we are obtaining rdf of same moleculetype, we have not yet
		!stored its COM coordinates. 
		!moleculetype2 COM coordinates is the same as moleculetype1 COM coordinates

		IF (species1_ndx == species2_ndx) THEN
			mtype2_rxcom = mtype1_rxcom
			mtype2_rycom = mtype1_rycom
			mtype2_rzcom = mtype1_rzcom
		END IF

		!still need to divide the com by the total mass
		mtype1_rxcom(:) = REAL(mtype1_rxcom(:)/mtype1_MW)
		mtype1_rycom(:) = REAL(mtype1_rycom(:)/mtype1_MW)
		mtype1_rzcom(:) = REAL(mtype1_rzcom(:)/mtype1_MW)
		mtype2_rxcom(:) = REAL(mtype2_rxcom(:)/mtype2_MW)
		mtype2_rycom(:) = REAL(mtype2_rycom(:)/mtype2_MW)
		mtype2_rzcom(:) = REAL(mtype2_rzcom(:)/mtype2_MW)
		
		!Start binning distances
		IF (species1_ndx == species2_ndx) THEN
			IF (i>=begin_frame) THEN
				DO m=0,nmolecules(i,species1_ndx)-2

					DO n=m+1,nmolecules(i,species1_ndx)-1

						rxijp = mtype1_rxcom(m)-mtype2_rxcom(n)  
						ryijp = mtype1_rycom(m)-mtype2_rycom(n)  
						rzijp = mtype1_rzcom(m)-mtype2_rzcom(n)  
						!obtain minimum image distances
						CALL MIN_IMAGE(Lx(i),Ly(i),Lz(i),rxijp,ryijp,rzijp,rxij,ryij,rzij)
						rsq = SQRT(rxij*rxij+ryij*ryij+rzij*rzij)
						IF (rsq > Lmin/2.0) THEN
							CONTINUE
						ELSE
							bin_ndx = NINT(rsq/binwidth)
							ncounts(bin_ndx) = ncounts(bin_ndx) + 1
						ENDIF
						!temporarily store counts  of distances that lie within each bin
					END DO
				END DO

				!average counts and store in rdf
				!We also need to make sure that the counts are normalized by the
				!number of particles
				DO p=0,nbins-1
					rdf(p) = rdf(p) + 2.0*REAL(ncounts(p))/REAL(nmolecules(i,species1_ndx))
				END DO
				!reset number of counts for next frame
				ncounts(:) =0
				!obtain bulk density of molecules in box
				bulk_density = bulk_density + &
							nmolecules(i,species1_ndx)/(Lx(i)*Ly(i)*Lz(i))
			END IF
		ELSE
			IF (i>=begin_frame) THEN
				DO m=0,nmolecules(i,species1_ndx)-1

					DO n=0,nmolecules(i,species2_ndx)-1
						rxijp = mtype1_rxcom(m)-mtype2_rxcom(n)
						ryijp = mtype1_rycom(m)-mtype2_rycom(n)
						rzijp = mtype1_rzcom(m)-mtype2_rzcom(n)
						!obtain minimum image distances
						CALL MIN_IMAGE(Lx(i),Ly(i),Lz(i),rxijp,ryijp,rzijp,rxij,ryij,rzij)
						rsq = SQRT(rxij*rxij+ryij*ryij+rzij*rzij)
						!print *, rsq
						IF (rsq > Lmin/2.0) THEN
							CONTINUE
						ELSE
							bin_ndx = NINT(rsq/binwidth)
							ncounts(bin_ndx) = ncounts(bin_ndx) + 1
						ENDIF
	
					END DO
				END DO
				!average counts and store in rdf
				!We also need to make sure that the counts are normalized by the
				!number of particles
				DO p=0,nbins-1
					rdf(p) = rdf(p) &
					+ REAL(ncounts(p))
				END DO
				!reset number of counts for next frame
				ncounts(:) =0
				!obtain bulk density of molecules in box
				bulk_density = bulk_density + &
							(nmolecules(i,species1_ndx) &
							* nmolecules(i,species2_ndx))/(Lx(i)*Ly(i)*Lz(i))
			END IF
		END IF		


		
	END DO
	!divide by number of frames
	rdf = rdf/(end_frame-begin_frame)
	bulk_density = bulk_density/(end_frame-begin_frame)



!CALCULATE ATOM-ATOM RDF
ELSE
	ALLOCATE(atype1_rx(0:max_nmols-1), atype1_ry(0:max_nmols-1), atype1_rz(0:max_nmols-1), &
		atype2_rx(0:max_nmols-1), atype2_ry(0:max_nmols-1), atype2_rz(0:max_nmols-1))
	!read xyz file lines for each frame
	DO i=0,end_frame-1
		!reset counter for atype index
		ia = 0
		ja = 0
		!first 2 lines is molecule #, and newline
		READ(*,*)
		READ(*,*)
		!loop over number of species
		DO j=0,nspecies-1
			!loop over number of molecules to store coordinates of atomtype
			DO k=0,nmolecules(i,j)-1
				DO l=0,natoms(j)-1
					READ (*,*) atype, rxp,ryp,rzp
					IF (atype1_name == atype .AND. atype1_ndx == l) THEN
						atype1_rx(ia) = rxp
						atype1_ry(ia) = ryp
						atype1_rz(ia) = rzp
						ia = ia+1
					ELSEIF (atype2_name == atype .AND. atype2_ndx == l) THEN
						atype2_rx(ja) = rxp
						atype2_ry(ja) = ryp
						atype2_rz(ja) = rzp
						ja= ja+1
					END IF
				END DO
			END DO
		
		END DO

		!if we are obtaining rdf of same atomtype, we have not yet
		!store its coordinates. 
		!atomtype 2 coordinates is the same as atomtype1 coordinates
		IF (atype1_name == atype2_name .AND. atype1_ndx == atype2_ndx) THEN
			atype2_rx = atype1_rx
			atype2_ry = atype1_ry
			atype2_rz = atype1_rz
			ja = ia
		END IF

   	!Start binning distances
		IF (species1_ndx == species2_ndx) THEN
			IF (i>=begin_frame) THEN
				DO m=0,ia-2

					DO n=m+1,ja-1
						rxijp = atype1_rx(m)-atype2_rx(n)
						ryijp = atype1_ry(m)-atype2_ry(n)
						rzijp = atype1_rz(m)-atype2_rz(n)
						!obtain minimum image distances
						CALL MIN_IMAGE(Lx(i),Ly(i),Lz(i),rxijp,ryijp,rzijp,rxij,ryij,rzij)
						rsq = SQRT(rxij*rxij+ryij*ryij+rzij*rzij)
						IF (rsq > Lmin/2.0) THEN
							CONTINUE
						ELSE
							bin_ndx = NINT(rsq/binwidth)
							ncounts(bin_ndx) = ncounts(bin_ndx) + 1
						ENDIF
						!temporarily store counts  of distances that lie within each bin
					END DO
				END DO

				!average counts and store in rdf
				!We also need to make sure that the counts are normalized by the
				!number of particles
				DO p=0,nbins-1
					rdf(p) = rdf(p) + 2.0*REAL(ncounts(p))/REAL(nmolecules(i,species1_ndx))
				END DO
				!reset number of counts for next frame
				ncounts(:) =0
				!obtain bulk density of molecules in box
				bulk_density = bulk_density + &
							nmolecules(i,species1_ndx)/(Lx(i)*Ly(i)*Lz(i))

			END IF
		ELSE
			IF (i>=begin_frame) THEN
				DO m=0,ia-1

					DO n=0,ja-1
						rxijp = atype1_rx(m)-atype2_rx(n)
						ryijp = atype1_ry(m)-atype2_ry(n)
						rzijp = atype1_rz(m)-atype2_rz(n)
						!obtain minimum image distances
						CALL MIN_IMAGE(Lx(i),Ly(i),Lz(i),rxijp,ryijp,rzijp,rxij,ryij,rzij)
						rsq = SQRT(rxij*rxij+ryij*ryij+rzij*rzij)
						IF (rsq > Lx(i)/2.0) THEN
							CONTINUE
						ELSE
							bin_ndx = NINT(rsq/binwidth)
							ncounts(bin_ndx) = ncounts(bin_ndx) + 1
						ENDIF
	
					END DO
				END DO
				!average counts and store in rdf
				!We also need to make sure that the counts are normalized by the
				!number of particles
				DO p=0,nbins-1
					rdf(p) = rdf(p) &
					+ REAL(ncounts(p))
				END DO
				!reset number of counts for next frame
				ncounts(:) =0
				!obtain bulk density of molecules in box
				bulk_density = bulk_density + &
							(nmolecules(i,species1_ndx) &
							* nmolecules(i,species2_ndx))/(Lx(i)*Ly(i)*Lz(i))
			END IF
		END IF		


   	
	END DO
   !divide by number of frames
	rdf = rdf/(end_frame-begin_frame)
	bulk_density = bulk_density/(end_frame-begin_frame)
	!we still need to normalize the rdf by the shell volumes
END IF

CLOSE(unit = 5)


shell_volume = 0.0
DO i=0,nbins-1
	rdf(i) = rdf(i)/(((shell_volume+binwidth)**3 - shell_volume**3) *4.0/3.0 * 3.1415*bulk_density)
	shell_volume = shell_volume+binwidth
END DO


END SUBROUTINE

!SUBROUTINE CAS_ANGLE(xyzfile,nspecies, &
!					nframes,begin_frame,end_frame, &
!					Lx,Ly,Lz,Lmin, &
!					species1_ndx, species2_ndx, species3_ndx, &
!					atype1_name, atype2_name, atype3_name, &
!					atype_mass,max_natoms,max_nmols,natoms, &
!					nmolecules,com_flag, angle)
!!**************************************************************************
!!Routine to generate the radial distribution function for given pair indices
!!
!! 		Input Variables: 
!!		xyzfile (INT) - coordinate file containing all frames
!!		nbins (INT) - number of bins
!!		binwidth (REAL) - based on minimum box length
!!		nspecies (INT) - number of species contained in the system
!!		nframes (INT) - number of frames we are analyzing
!!		begin_frame (INT) - beginning frame index
!!		end_frame (INT) - ending_frame index
!!		Lx, Ly, Lz (ARRAY) - box lengths 
!!		Lmin (REAL) - global minimum box length from Lx,Ly, and Lz
!!		species1_ndx, species2_ndx (INT) - species index
!!		atype1_name, atype_name (CHAR) -name of atomtype
!!		atype_mass (MATRIX) - contains masses of each atomtype for each species
!!		max_natoms (INT) - max number of atoms out of atoms in species
!!		max_nmols (INT) - max number of molecules out of molecules in species
!!		natoms (ARRAY) - number of atoms in each species
!!		nmolecules (MATRIX) - array containing number of molecules for each 
!!		species in each frame
!!		com_flag (LOGICAL) - flag determining analysis: com-com RDF or atom-atom RDF
!!		angle (output) - rdf array
!!
!!**************************************************************************
!IMPLICIT NONE
!CHARACTER(120), INTENT(IN) :: xyzfile
!INTEGER, INTENT(IN) :: nslices, nspecies, nframes, max_natoms, begin_frame, end_frame
!REAL, INTENT(IN), DIMENSION(0:nframes-1) :: Lx, Ly, Lz
!!f2py depend(nframes) :: Lx, Ly, Lz
!REAL,INTENT(IN), DIMENSION(0:max_natoms-1,0:nspecies-1) :: atype_mass
!!f2py depend(max_natoms,nspecies) :: atype_mass
!INTEGER, INTENT(IN), DIMENSION(0:nspecies-1) :: natoms
!!f2py depend(nspecies) :: natoms
!INTEGER, INTENT(IN), DIMENSION(0:nframes-1,0:nspecies-1) :: nmolecules
!!f2py depend(nframes,nspecies) :: nmolecules
!REAL, INTENT(INOUT), DIMENSION(0:nslices-1) :: angle
!!f2py intent(in,out) :: angle
!!f2py depend(nslices) :: angle
!
!INTEGER :: i,j,k,l, slice
!REAL :: rxp,ryp,rzp, binwidth, binvol,conv
!CHARACTER(6) :: atype
!
!conv = 1660.54!amu/Angstrom^3 to kg/m^3 conversion factor
!
!OPEN(unit=5,file=xyzfile)
!
!!loop over #frames
!DO i=0,end_frame-1
!	!first 2 lines is molecule #, and newline
!	READ(*,*)
!	READ(*,*)
!
!	!define binwidth in the frame
!	binwidth = REAL(Lz(i)/nslices)
!	binvol = REAL(nslices/(Lx(i)*Ly(i)*Lz(i)))
!	!loop over number of species
!	DO j=0,nspecies-1
!		!loop over number of molecules
!		DO k=0,nmolecules(i,j)-1
!			DO l=0, natoms(j)-1 
!				READ (*,*) atype, rxp,ryp,rzp
!				IF (i>=begin_frame) THEN
!					slice = ABS(NINT((rzp+Lz(i)/2.0)/binwidth))
!					density(slice) = density(slice)+1.0*binvol *atype_mass(l,j)
!				END IF
!			END DO
!		END DO
!
!	END DO
!
!END DO
!
!CLOSE(unit = 5)
!angle = angle/(end_frame-begin_frame)*conv
!
!END SUBROUTINE




SUBROUTINE MIN_IMAGE(Lx,Ly,Lz,rxijp,ryijp,rzijp,rxij,ryij,rzij)
!**************************************************************************
!Routine to obtain the minimum image distance
!       Input Variables: 
!
!		Lx, Ly, Lz (REAL) - instantaneous box lengths
!		rxijp, ryijp, rzijp (REAL) - pairwise distances
!		rxij, ryij, rzij (REAL) - minimum image pairwise distance
!**************************************************************************
REAL, INTENT(IN) :: Lx,Ly,Lz, rxijp, ryijp, rzijp 
REAL, INTENT(OUT) :: rxij, ryij, rzij
REAL :: hLx, hLy, hLz


hLx = Lx/2.0
hLy = Ly/2.0
hLz = Lz/2.0

rxij = rxijp
ryij = ryijp
rzij = rzijp

IF (rxijp > hLx) THEN
	rxij = rxijp - Lx
ELSEIF (rxijp < -hLx) THEN
	rxij = rxijp + Lx
ENDIF


IF (ryijp > hLy) THEN
	ryij = ryijp - Ly
ELSEIF (ryijp < -hLy) THEN
	ryij = ryijp + Ly
ENDIF


IF (rzijp > hLz) THEN
	rzij = rzijp - Lz
ELSEIF (rzijp < -hLz) THEN
	rzij = rzijp + Lz
ENDIF

END SUBROUTINE





