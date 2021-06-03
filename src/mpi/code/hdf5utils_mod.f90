! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://www.hdfgroup.org/licenses.               *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! Fortran parallel example.  Copied from Tutorial's example program of
! dataset.f90.
module hdf5utils_mod

    use, intrinsic :: iso_fortran_env
    implicit none
    ! integer, parameter :: dp = real64

    contains
     subroutine writehdf5(rays)
        USE HDF5 ! This module contains all necessary modules 
        USE MPI

        IMPLICIT NONE
        LOGICAL :: debug = .false.
        INTEGER :: rays
        CHARACTER(LEN=10), PARAMETER :: filename = "sds_col.h5"  ! File name
        CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

        INTEGER(HID_T) :: file_id       ! File identifier 
        INTEGER(HID_T) :: dset_id       ! Dataset identifier 
        INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
        INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
        INTEGER(HID_T) :: plist_id      ! Property list identifier 

        INTEGER(HSIZE_T), DIMENSION(2) :: dimsf ! Dataset dimensions.
    !     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi

        INTEGER(HSIZE_T), DIMENSION(2) :: count  
        INTEGER(HSSIZE_T), DIMENSION(2) :: offset 
        REAL(kind(0.d0)), ALLOCATABLE :: data (:,:)  ! Data to write
        INTEGER :: rank = 2 ! Dataset rank 

        INTEGER :: error, error_n  ! Error flags
        !
        ! MPI definitions and calls.
        !
        INTEGER :: mpierror       ! MPI error flag
        INTEGER :: comm, info
        INTEGER :: mpi_size, mpi_rank

        dimsf = (/7,rays/)
        dimsfi = (/7,rays/)

        comm = MPI_COMM_WORLD
        info = MPI_INFO_NULL
        CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
        CALL MPI_COMM_RANK(comm, mpi_rank, mpierror) 
        !
        ! Initialize FORTRAN predefined datatypes
        !
        if(debug) write(6,*) "Rank ",mpi_rank," opens Fortran interface..."
        CALL h5open_f(error) 

        ! 
        ! Setup file access property list with parallel I/O access.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," creates P I/O access..."

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

        !
        ! Create the file collectively.
        ! 
        if(debug) write(6,*) "Rank ",mpi_rank," creates file collectively..."

        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
        CALL h5pclose_f(plist_id, error)
        !
        ! Create the data space for the  dataset. 
        !
        if(debug) write(6,*) "Rank ",mpi_rank," creates data space..."

        CALL h5screate_simple_f(rank, dimsf, filespace, error)

        !
        ! Create the dataset with default properties.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," creates dataset properties..."

        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
                        dset_id, error)
        CALL h5sclose_f(filespace, error)
        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file. 
        !
        count(1) = dimsf(1)
        count(2) = dimsf(2)/mpi_size 
        offset(1) = 0
        offset(2) = mpi_rank * count(2) 
        CALL h5screate_simple_f(rank, count, memspace, error) 
        ! 
        ! Select hyperslab in the file.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," creates hyperslab..."

        CALL h5dget_space_f(dset_id, filespace, error)
        CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)
        ! 
        ! Initialize data buffer with trivial data.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," creates data..."

        ALLOCATE ( data(count(1),count(2)))
        data = mpi_rank + 10
        !
        ! Create property list for collective dataset write
        !
        if(debug) write(6,*) "Rank ",mpi_rank," creates properties for data write..."

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        
        if(debug) write(6,*) "Rank ",mpi_rank," writes file..."

        !
        ! Write the dataset collectively. 
        !
        ! CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsfi, error, &
        !                 file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
        !
        ! Write the dataset independently. 
        !

        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsfi, error, &
                        file_space_id = filespace, mem_space_id = memspace)
        !
        ! Deallocate data buffer.
        !
                        if(debug) write(6,*) "Rank ",mpi_rank," deallocates data..."

        DEALLOCATE(data)
        !
        ! Close dataspaces.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," closes dataspaces..."

        CALL h5sclose_f(filespace, error)
        CALL h5sclose_f(memspace, error)

        !
        ! Close the dataset and property list.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," closes property lists..."

        CALL h5dclose_f(dset_id, error)
        CALL h5pclose_f(plist_id, error)

        !
        ! Close the file.
        !
        if(debug) write(6,*) "Rank ",mpi_rank," closes file handle..."

        CALL h5fclose_f(file_id, error)

        !
        ! Close FORTRAN predefined datatypes.
        !
        CALL h5close_f(error)
        if(debug) write(6,*) "Rank ",mpi_rank," finishes writes..."

    
    end subroutine 

end module

! USE HDF5 ! This module contains all necessary modules
!         USE MPI
!         IMPLICIT NONE
!         INTEGER :: maxrayp
!         CHARACTER(LEN=10), PARAMETER :: default_fname = "sds.h5"  ! Default name
!         CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

!         CHARACTER(LEN=100) :: filename  ! File name
!         INTEGER        :: fnamelen	     ! File name length
!         INTEGER(HID_T) :: file_id       ! File identifier
!         INTEGER(HID_T) :: dset_id       ! Dataset identifier
!         INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
!         INTEGER(HID_T) :: plist_id      ! Property list identifier

!         INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/5,8/) ! Dataset dimensions.
!     !     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
!     !     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/5,8/)
!         INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi

!         REAL(kind(0.d0)), ALLOCATABLE :: data(:,:)   ! Data to write
!         INTEGER :: rank = 2 ! Dataset rank

!         INTEGER :: error, error_n  ! Error flags
!         INTEGER :: i, j
!         !
!         ! MPI definitions and calls.
!         !
!         INTEGER :: mpierror       ! MPI error flag
!         INTEGER :: comm, info
!         INTEGER :: mpi_size, mpi_rank
!         comm = MPI_COMM_WORLD
!         info = MPI_INFO_NULL
!         ! CALL MPI_INIT(mpierror)
!         CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
!         CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)
!         !
!         ! Initialize data buffer with trivial data.
!         !
!         ALLOCATE ( data(dimsf(1),dimsf(2)))
!         do i = 1, dimsf(2)
!         do j = 1, dimsf(1)
!             data(j,i) = j**2 - 1 + (i-1)*dimsf(1)**2
!         enddo
!         enddo
!         !
!         ! Initialize FORTRAN interface
!         !
!         CALL h5open_f(error)

!         !
!         ! Setup file access property list with parallel I/O access.
!         !
!         CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
!         CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

!         !
!         ! Figure out the filename to use.  If your system does not support
!         ! getenv, comment that statement with this,
!         ! filename = ""
!         CALL get_environment_variable("HDF5_PARAPREFIX", filename)
!         fnamelen = LEN_TRIM(filename)
!         if ( fnamelen == 0 ) then
!         filename = default_fname
!         else
!         filename = filename(1:fnamelen) // "/" // default_fname
!         endif
!         print *, "Using filename = ", filename

!         !
!         ! Create the file collectively.
!         !
!         CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
!         CALL h5pclose_f(plist_id, error)
!         !
!         ! Create the data space for the  dataset.
!         !
!         CALL h5screate_simple_f(rank, dimsf, filespace, error)

!         !
!         ! Create the dataset with default properties.
!         !
!         CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
!                         dset_id, error)
!         !
!         ! Create property list for collective dataset write
!         !
!         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!         ! CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!         !
!         ! For independent write use
!         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
!         !

!         !
!         ! Write the dataset collectively.
!         !
!         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsfi, error, &
!                         xfer_prp = plist_id)
!         !
!         ! Deallocate data buffer.
!         !
!         DEALLOCATE(data)

!         !
!         ! Close resources.
!         !
!         CALL h5sclose_f(filespace, error)
!         CALL h5dclose_f(dset_id, error)
!         CALL h5pclose_f(plist_id, error)
!         CALL h5fclose_f(file_id, error)
!         ! Attempt to remove the data file.  Remove the line if the compiler
!         ! does not support it.
!         !CALL unlink(filename)

!         !
!         ! Close FORTRAN interface
!         !
!         CALL h5close_f(error)
                