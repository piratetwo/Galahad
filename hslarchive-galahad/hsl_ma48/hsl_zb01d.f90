! COPYRIGHT (c) 2007 Science & Technology Facilities Council
! Original date 1 August 2007. Version 1.0.0. (Sue Dollar)
!
! 1 Decemeber 2010 Version 2.0.0 (Jonathan Hogg)
!    Modify interface to allow specificaiton of source and destination ranges,
!    substantially rewrite and simplify code, add long integer support, remove
!    support for unallocated arrays on input

! To convert to single:
!    s/double/single
!    change myreal definition
! To convert to integer:
!    s/double/integer
!    s/real(myreal)/integer(myinteger)
!    s/REAL (kind=myreal)/INTEGER (kind=myinteger)
!    change myreal definintion to myinteger
! To convert to long integer:
!    s/double/long_integer
!    s/real(myreal)/integer(myinteger)
!    s/REAL (kind=myreal)/INTEGER (kind=myinteger)
!    change myreal definintion to myinteger
! To convert to double complex:
!    s/double/double_complex
!    s/COMPLEX (kind=mycomplex)/COMPLEX (kind=mycomplex)
!    s/complex(mycomplex)/complex(mycomplex)
!    change myreal definition to mycomplex
! To convert to complex:
!    s/double/complex
!    s/COMPLEX (kind=mycomplex)/COMPLEX (kind=mycomplex)
!    s/complex(mycomplex)/complex(mycomplex)
!    change myreal definition to mycomplex
MODULE hsl_zb01_double

   IMPLICIT NONE
   PRIVATE

   ! ---------------------------------------------------
   ! Precision
   ! ---------------------------------------------------

   INTEGER, PARAMETER :: myreal = kind(1.0d0)
   INTEGER, PARAMETER :: myint = kind(1)
   INTEGER, PARAMETER :: long = selected_int_kind(18)

   ! ---------------------------------------------------
   ! Error flags
   ! ---------------------------------------------------
   INTEGER (kind=myint), PARAMETER :: &
      zb01_err_lw = -1, &        ! lw<=lkeep/size(w) on input
      zb01_err_lw_1 = -2, &      ! lw<1 on input
      zb01_err_lw_both = -3, &   ! both -1 and -2
      zb01_err_lkeep = -4, &     ! lkeep>size(w) on input
      zb01_err_lkeep_l = -5, &   ! lkeep<1
      zb01_err_lkeep_both = -6, &! both -4 and -5
      zb01_err_filename = -7, &  ! filename too long
      zb01_err_file_size = -8, & ! file_size <2**12
      zb01_err_filename_exists = -9, & ! filename already exists
      zb01_err_mode = -10, &     ! mode out of range
      zb01_err_memory_alloc = -11, & ! memory alloc error
      zb01_err_memory_dealloc = -12, & ! memory dealloc error
      zb01_err_inquire = -13, &  ! error in Fortran inquire statement
      zb01_err_open = -14, &     ! error in Fortran open statement
      zb01_err_read = -15, &     ! error in Fortran read statement
      zb01_err_write = -16, &    ! error in Fortran write statement
      zb01_err_close = -17, &    ! error in Fortran close statement
      zb01_err_src_dest = -18, & ! src and dest sizes do not match
      zb01_err_w_unalloc = -19   ! w is unallocated on entry

   ! ---------------------------------------------------
   ! Warning flags
   ! ---------------------------------------------------
   INTEGER (kind=myint), PARAMETER :: &
      zb01_warn_lw = 1 ! size_out changed

   ! ---------------------------------------------------
   ! Derived type definitions
   ! ---------------------------------------------------

   TYPE, PUBLIC :: zb01_info
      INTEGER :: flag = 0        ! error/warning flag
      INTEGER :: iostat = 0      ! holds Fortran iostat parameter
      INTEGER :: stat = 0        ! holds Fortran stat parameter
      INTEGER :: files_used = 0  ! unit number scratch file written to
   END TYPE zb01_info

   INTERFACE zb01_resize1
      MODULE PROCEDURE zb01_resize1_double
   END INTERFACE

   INTERFACE zb01_resize2
      MODULE PROCEDURE zb01_resize2_double
   END INTERFACE

   PUBLIC zb01_resize1, zb01_resize2

CONTAINS

SUBROUTINE zb01_resize1_double(w,size_in,size_out,info,src,dest,filename, &
      file_size,mode)

  ! --------------------------------------
  ! Expands w to have larger size and copies all or part of the
  ! original array into the new array
  ! --------------------------------------

  ! w: is a REAL allocatable array of INTENT(INOUT). Must be allocated on entry.
  REAL (kind=myreal), DIMENSION (:), ALLOCATABLE, INTENT (INOUT) :: w

  ! size_in: is an INTEGER OF INTENT(IN). It holds the extent of w on entry.
  INTEGER (kind=long), INTENT(IN) :: size_in

  ! size_out: is an INTEGER of INTENT(INOUT). It holds the required size of w
  ! on input and the actual size of w on output
  INTEGER (kind=long), INTENT (INOUT) :: size_out

  ! info: is of derived type ZB01_info with intent(out).
  ! info%flag = ZB01_ERR_LW if 1<=size_out<=lkeep/size(w) on input
  ! ZB01_ERR_LW_L if size_out<=0
  ! ZB01_ERR_LKEEP if lkeep>size(w) on input
  ! ZB01_ERR_LKEEP_L if lkeep<1 on input
  ! ZB01_ERR_FILENAME if filename too long
  ! ZB01_ERR_FILE_SIZE if file_size <2**12
  ! ZB01_ERR_FILENAME_EXISTS if filename already exists
  ! ZB01_ERR_MODE if mode out of range
  ! ZB01_ERR_MEMORY_ALLOC if memory alloc error
  ! ZB01_ERR_MEMORY_DEALLOC if memory dealloc error
  ! ZB01_ERR_INQUIRE if error in Fortran inquire statement
  ! ZB01_ERR_OPEN if error in Fortran open statement
  ! ZB01_ERR_READ if error in Fortran read statement
  ! ZB01_ERR_WRITE if error in Fortran write statement
  ! ZB01_ERR_CLOSE if error in Fortran close statement
  ! ZB01_WARN_LW if size_out changed by subroutine
  ! ZB01_WARN_W if w not allocated on input
  ! info%iostat holds Fortran iostat parameter
  ! info%stat holds Fortran stat parameter
  ! info%unit holds unit number scratch file written to (negative if
  ! not)
  TYPE (zb01_info), INTENT (OUT) :: info

  ! src and dest: are OPTIONAL INTEGER arrays of INTENT(IN). They specify
  ! the source and destination ranges for any values to be kept.
  ! dest(2)-dest(1) must equal src(2)-src(1) and minval(src,dest)>0.
  INTEGER (kind=long), DIMENSION(2), INTENT (IN), OPTIONAL :: src
  INTEGER (kind=long), DIMENSION(2), INTENT (IN), OPTIONAL :: dest

  ! filename: is an OPTIONAL STRING OF CHARACTERS of INTENT(IN).
  ! It holds the name of the file to be used
  CHARACTER (len=*), INTENT (IN), OPTIONAL :: filename

  ! file_size: is an OPTIONAL INTEGER of INTENT(IN). It holds the length
  ! of
  ! of the files to be used if necessary
  INTEGER (kind=long), INTENT (IN), OPTIONAL :: file_size

  ! mode: is an OPTIONAL INTEGER of INTENT(IN). If mode==0, then the
  ! subroutine will firstly try to use a temporary array. If mode==1,
  ! then
  ! the subroutine will immediately use temporary files
  INTEGER, INTENT (IN), OPTIONAL :: mode

  ! --------------------------------------
  ! Local variables
  ! --------------------------------------

  ! wtemp: is a REAL allocatable array
  REAL (kind=myreal), DIMENSION (:), ALLOCATABLE :: wtemp

  ! length to use for wtemp
  INTEGER :: lwtemp

  ! Fortran stat parameter
  INTEGER :: stat

  INTEGER :: number_files, mode_copy
  INTEGER (kind=long) :: file_size_copy
  INTEGER, DIMENSION (:), ALLOCATABLE :: units
  integer (kind=long), dimension(2,2) :: src_copy, dest_copy

  ! --------------------------------------
  ! Check input for errors
  ! --------------------------------------

  info%flag = 0

  ! Check w is allocated
  if(.not.allocated(w)) then
    info%flag = zb01_err_w_unalloc
    return
  endif

  ! Setup src_copy and dest_copy
  src_copy(1,1) = 1
  src_copy(2,1) = size_in
  src_copy(1,2) = 1
  src_copy(2,2) = 1
  if(present(src)) then
    src_copy(:,1) = src(:)
  elseif(present(dest)) then
    src_copy(:,1) = dest(:)
  endif
  dest_copy(:,:) = src_copy(:,:)
  if(present(dest)) then
    dest_copy(:,1) = dest(:)
  endif

  if(minval(src_copy).lt.0 .or. minval(dest_copy).lt.0) &
    info%flag = zb01_err_lkeep_l

  if(maxval(src_copy) .gt. size_in) then
    if(info%flag.eq.zb01_err_lkeep_l) then
      info%flag = zb01_err_lkeep_both
    else
      info%flag = zb01_err_lkeep
    endif
  endif
  if(info%flag.lt.0) return
  
  if(any(src_copy(2,:)-src_copy(1,:) .ne. &
      dest_copy(2,:)-dest_copy(1,:))) then
    info%flag = zb01_err_src_dest
    return
  endif

  if(size_out.lt.1) then
    info%flag = zb01_err_lw_1
    return
  endif
  if(maxval(dest_copy).gt.size_out) then
    info%flag = zb01_err_lw
    return
  endif

  IF (present(filename)) THEN
    IF (len(filename)>400) THEN
      ! Error: len(filename) > 400
      info%flag = zb01_err_filename
      RETURN
    END IF
  END IF


  file_size_copy = 2**22_long
  if (present(file_size)) file_size_copy = file_size

  IF (file_size_copy<2**12) THEN
    ! Error: filesize < 2**12
    info%flag = zb01_err_file_size
    RETURN
  END IF

  mode_copy = 0
  if(present(mode)) mode_copy = mode
  IF (mode_copy<0 .OR. mode_copy>1) THEN
    ! Error: filesize < 2**12
    info%flag = zb01_err_mode
    RETURN
  END IF

  IF (src_copy(2,1)-src_copy(1,1)+1.le.0) THEN

    ! --------------------------------------
    ! No entries need to be copied
    ! --------------------------------------

    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out==size_in) RETURN

    ! --------------------------------------
    ! Deallocate w
    ! --------------------------------------
    DEALLOCATE (w,STAT=info%stat)
    IF (info%stat>0) THEN
      info%flag = zb01_err_memory_dealloc
      RETURN
    END IF

    ! --------------------------------------
    ! Reallocate w
    ! --------------------------------------
    ALLOCATE (w(size_out),STAT=info%stat)
    IF (info%stat>0) THEN
      ! --------------------------------------
      ! Allocation of w failed
      ! --------------------------------------
      info%flag = zb01_err_memory_alloc
      RETURN
    END IF

  ELSE

    ! --------------------------------------
    ! entries need to be copied
    ! --------------------------------------


    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out==size_in .and. all(src_copy(:,:).eq.dest_copy(:,:))) RETURN

    ! --------------------------------------
    ! Attempt to allocate temporary array
    ! --------------------------------------
    stat = 0
    ! Length of temporary array
    lwtemp = src_copy(2,1) - src_copy(1,1) + 1

    ! Allocate temporary array
    IF (mode_copy==0) ALLOCATE (wtemp(lwtemp),STAT=stat)

    IF (stat==0 .AND. mode_copy==0) THEN
      ! --------------------------------------
      ! Allocation successful
      ! --------------------------------------

      ! --------------------------------------
      ! Copy entries into wtemp
      ! --------------------------------------
      wtemp(1:lwtemp) = w(src_copy(1,1):src_copy(2,1))
      !print *, "copy w to wtemp"

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out),STAT=stat)
      IF (stat>0) THEN
        !print *, "fail alloc of w, copy wtemp to file", src_copy, size_out


        ! --------------------------------------
        ! Allocation not successful
        ! --------------------------------------
        call write_to_file(wtemp, size_in, src_copy, &
          units, number_files, file_size_copy, info, filename)
        if(info%flag.lt.0) return

        ! --------------------------------------
        ! Deallocate wtemp
        ! --------------------------------------
        DEALLOCATE (wtemp,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = zb01_err_memory_dealloc
          RETURN
        END IF

        ! --------------------------------------
        ! Reallocate w
        ! --------------------------------------
        ALLOCATE (w(size_out),STAT=info%stat)
        IF (info%stat>0) THEN
          !print *, "still can't alloc w. give up (try and preserve)."
          ! --------------------------------------
          ! Allocation of w to desired size failed
          ! Try size=lwtemp
          ! --------------------------------------


          ALLOCATE (w(lwtemp),STAT=info%stat)
          IF (info%stat==0) THEN
            ! --------------------------------------
            ! Reassign lw
            ! --------------------------------------
            size_out = lwtemp
            info%flag = zb01_warn_lw

          ELSE
            !print *, "can't even do that!"

            ! --------------------------------------
            ! Allocation of w still failed - delete files
            ! --------------------------------------
            call delete_files(units, number_files, info, filename)
            if(info%flag.lt.0) return

            info%flag = zb01_err_memory_alloc
            RETURN
          END IF
        END IF

        !print *, "restore entries to w"
        ! --------------------------------------
        ! Copy entries back into w
        ! --------------------------------------
        call read_from_file(w, size_in, dest_copy, units, &
          number_files, file_size_copy, info, filename)

        RETURN
      END IF

      !print *, "copy wtemp back to w"
      ! --------------------------------------
      ! Copy entries
      ! --------------------------------------
      w(dest_copy(1,1):dest_copy(2,1)) = wtemp(1:lwtemp)

      ! --------------------------------------
      ! Deallocate wtemp
      ! --------------------------------------
      DEALLOCATE (wtemp,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF
      RETURN

    ELSE
      !print *, "can't alloc wtemp, stow direct in file"
      ! --------------------------------------
      ! Allocation not successful
      ! --------------------------------------
      call write_to_file(w, size_in, src_copy, units, &
        number_files, file_size_copy, info, filename)
      if(info%flag.lt.0) return

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out),STAT=info%stat)
      IF (info%stat>0) THEN
        !print *, "can't alloc w"
        ! --------------------------------------
        ! Allocation of w to desired size failed
        ! Try size=lwtemp
        ! --------------------------------------

        ALLOCATE (w(lwtemp),STAT=info%stat)
        IF (info%stat==0) THEN
          ! --------------------------------------
          ! Reassign lw
          ! --------------------------------------
          size_out = lwtemp
          info%flag = zb01_warn_lw

        ELSE
          ! --------------------------------------
          ! Allocation of w failed again - delete files
          ! --------------------------------------
          call delete_files(units, number_files, info, filename)
          if(info%flag.lt.0) return

          info%flag = zb01_err_memory_alloc
          RETURN
        END IF
      END IF

      !print *, "restore values from file"
      call read_from_file(w, size_in, dest_copy, units, &
        number_files, file_size_copy, info, filename)

    END IF
  END IF

END SUBROUTINE zb01_resize1_double

! ---------------------------------------------------------------

SUBROUTINE zb01_resize2_double(w,size_in,size_out,info,src,dest,filename, &
      file_size,mode)

  ! --------------------------------------
  ! Expands w to have larger size and copies all or part of the
  ! original array into the new array
  ! --------------------------------------

  ! w: is a REAL allocatable array of INTENT(INOUT). Must be allocated on entry.
  REAL (kind=myreal), DIMENSION (:,:), ALLOCATABLE, INTENT (INOUT) :: w

  ! size_in: is an INTEGER array INTENT(IN). It holds the extent of w on entry.
  INTEGER (kind=long), INTENT(IN) :: size_in(2)

  ! size_out: is an INTEGER array of rank-one with size 2 and of
  ! INTENT(INOUT).
  ! On input, size_out(1) holds the required number of rows of w and size_out(2)
  ! holds
  ! the required number of columns of w. On successful output it
  ! contains
  ! the number of rows and columns that w has.
  INTEGER (kind=long), INTENT (INOUT) :: size_out(2)

  ! info: is of derived type ZB01_info with intent(out).
  ! info%flag = ZB01_ERR_LW if 1<=size_out<=lkeep/size(w) on input
  ! ZB01_ERR_LW_L if size_out<=0
  ! ZB01_ERR_LKEEP if lkeep>size(w) on input
  ! ZB01_ERR_LKEEP_L if lkeep<1 on input
  ! ZB01_ERR_FILENAME if filename too long
  ! ZB01_ERR_FILE_SIZE if file_size <2**12
  ! ZB01_ERR_FILENAME_EXISTS if filename already exists
  ! ZB01_ERR_MODE if mode out of range
  ! ZB01_ERR_MEMORY_ALLOC if memory alloc error
  ! ZB01_ERR_MEMORY_DEALLOC if memory dealloc error
  ! ZB01_ERR_INQUIRE if error in Fortran inquire statement
  ! ZB01_ERR_OPEN if error in Fortran open statement
  ! ZB01_ERR_READ if error in Fortran read statement
  ! ZB01_ERR_WRITE if error in Fortran write statement
  ! ZB01_ERR_CLOSE if error in Fortran close statement
  ! ZB01_WARN_LW if size_out changed by subroutine
  ! ZB01_WARN_W if w not allocated on input
  ! info%iostat holds Fortran iostat parameter
  ! info%stat holds Fortran stat parameter
  ! info%files_used holds unit number of files written to
  TYPE (zb01_info), INTENT (OUT) :: info

  ! src and dest: are OPTIONAL INTEGER arrays of INTENT(IN). They specify
  ! the source and destination ranges for any values to be kept.
  ! dest(2,:)-dest(1,:) must equal src(2,:)-src(1,:) and must be
  ! non-negative.
  INTEGER (kind=long), DIMENSION(2,2), INTENT (IN), OPTIONAL :: src
  INTEGER (kind=long), DIMENSION(2,2), INTENT (IN), OPTIONAL :: dest

  ! filename: is an OPTIONAL STRING OF CHARACTERS of INTENT(IN).
  ! It holds the name of the file to be used
  CHARACTER (len=*), INTENT (IN), OPTIONAL :: filename


  ! file_size: is an OPTIONAL INTEGER of INTENT(IN). It holds the length
  ! of the files to be used if necessary
  INTEGER (kind=long), INTENT (IN), OPTIONAL :: file_size

  ! mode: is an OPTIONAL INTEGER of INTENT(IN). If mode==0, then the
  ! subroutine will firstly try to use a temporary array. If mode==1,
  ! then the subroutine will immediately use temporary files
  INTEGER, INTENT (IN), OPTIONAL :: mode


  ! --------------------------------------
  ! Local variables
  ! --------------------------------------

  ! wtemp: is a REAL allocatable array
  REAL (kind=myreal), DIMENSION (:,:), ALLOCATABLE :: wtemp

  ! lengths to use for wtemp
  INTEGER :: lwtemp1, lwtemp2

  ! Fortran stat parameter
  INTEGER :: stat

  INTEGER :: number_files, mode_copy
  INTEGER (kind=long) :: file_size_copy
  INTEGER, DIMENSION (:), ALLOCATABLE :: units

  integer(long), dimension(2,2) :: src_copy, dest_copy

  info%flag = 0

  ! --------------------------------------
  ! Check input for errors
  ! --------------------------------------

  ! Check w is allocated
  if(.not.allocated(w)) then
    info%flag = zb01_err_w_unalloc
    return
  endif

  src_copy(1,1) = 1
  src_copy(2,1) = size_in(1)
  src_copy(1,2) = 1
  src_copy(2,2) = size_in(2)
  if(present(src)) then
    src_copy(:,:) = src(:,:)
  elseif(present(dest)) then
    src_copy(:,:) = dest(:,:)
  endif
  dest_copy(:,:) = src_copy(:,:)
  if(present(dest)) dest_copy(:,:) = dest(:,:)

  if(minval(src_copy).lt.0 .or. minval(dest_copy).lt.0) &
    info%flag = zb01_err_lkeep_l

  if(maxval(src_copy(:,1)).gt.size_in(1) .or. &
      maxval(src_copy(:,2)).gt.size_in(2)) then
    if(info%flag.eq.zb01_err_lkeep_l) then
      info%flag = zb01_err_lkeep_both
    else
      info%flag = zb01_err_lkeep
    endif
  endif
  if(info%flag.lt.0) return

  ! Note: be careful to only return one error for each dimension
  ! (combination of flags possible if dimensions are differently wrong)
  ! this is to match v1.0.0 behaviour
  if(minval(size_out).lt.1) info%flag = zb01_err_lw_1
  if((maxval(dest_copy(:,1)).gt.size_out(1) .and. size_out(1).ge.1) .or. &
     (maxval(dest_copy(:,2)).gt.size_out(2) .and. size_out(2).ge.1)) then
    if(info%flag.eq.zb01_err_lw_1) then
      info%flag = zb01_err_lw_both
    else
      info%flag = zb01_err_lw
    endif
  endif
  if(info%flag.lt.0) return
  
  if(any(src_copy(2,:)-src_copy(1,:) .ne. &
      dest_copy(2,:)-dest_copy(1,:))) then
    info%flag = zb01_err_src_dest
    return
  endif

  IF (present(filename)) THEN
    IF (len(filename)>400) THEN
      ! Error: len(filename) > 400
      info%flag = zb01_err_filename
      RETURN
    END IF
  END IF

  file_size_copy = 2**22_long
  if (present(file_size)) file_size_copy = file_size
  IF (file_size_copy<2**12) THEN
    ! Error: filesize < 2**12
    info%flag = zb01_err_file_size
    RETURN
  END IF

  mode_copy = 0
  if(present(mode)) mode_copy = mode
  IF (mode_copy<0 .OR. mode_copy>1) THEN
    ! Error: mode out of range
    info%flag = zb01_err_mode
    RETURN
  END IF

  IF (any(src_copy(2,:)-src_copy(1,:)+1.le.0)) THEN

    ! --------------------------------------
    ! No entries need to be copied
    ! --------------------------------------

    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out(1)==size_in(1) .AND. size_out(2)==size_in(2)) RETURN

    ! --------------------------------------
    ! Deallocate w
    ! --------------------------------------
    DEALLOCATE (w,STAT=info%stat)
    IF (info%stat>0) THEN
      info%flag = zb01_err_memory_dealloc
      RETURN
    END IF

    ! --------------------------------------
    ! Reallocate w
    ! --------------------------------------
    ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
    IF (info%stat>0) THEN
      ! --------------------------------------
      ! Allocation of w failed
      ! --------------------------------------
      info%flag = zb01_err_memory_alloc
      RETURN
    END IF

  ELSE

    ! --------------------------------------
    ! Entries need to be copied
    ! --------------------------------------


    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out(1)==size_in(1) .AND. size_out(2)==size_in(2) .and. &
         all(src_copy(:,:).eq.dest_copy(:,:))) RETURN

    ! --------------------------------------
    ! Attempt to allocate temporary array
    ! --------------------------------------
    stat = 0
    ! Size of temporary array
    lwtemp1 = src_copy(2,1) - src_copy(1,1) + 1
    lwtemp2 = src_copy(2,2) - src_copy(1,2) + 1

    ! Allocate temporary array
    IF (mode_copy==0) ALLOCATE (wtemp(lwtemp1,lwtemp2),STAT=stat)
    IF (stat==0 .AND. mode_copy==0) THEN
      ! --------------------------------------
      ! Allocation of temporary array successful
      ! --------------------------------------

      ! --------------------------------------
      ! Copy entries into wtemp
      ! --------------------------------------
      !print *, "copy to wtemp"
      wtemp(1:lwtemp1, 1:lwtemp2) = &
        w(src_copy(1,1):src_copy(2,1), src_copy(1,2):src_copy(2,2))

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
      IF (info%stat>0) THEN
        ! --------------------------------------
        ! Allocation of w not successful so copy temp into files
        ! --------------------------------------
        !print *, "copy wtemp to file"
        call write_to_file(wtemp, lwtemp1+0_long, src_copy, &
          units, number_files, file_size_copy, info, filename)
        if(info%flag.lt.0) return


        ! --------------------------------------
        ! Deallocate wtemp
        ! --------------------------------------
        DEALLOCATE (wtemp,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = zb01_err_memory_dealloc
          RETURN
        END IF

        ! --------------------------------------
        ! Reallocate w
        ! --------------------------------------
        ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
        IF (info%stat>0) THEN
          ! --------------------------------------
          ! Allocation of w to desired size failed
          ! Try size=lwtemp
          ! --------------------------------------

          ALLOCATE (w(lwtemp1,lwtemp2),STAT=info%stat)
          IF (info%stat==0) THEN
            ! --------------------------------------
            ! Reassign size_out
            ! --------------------------------------
            size_out(1) = lwtemp1
            size_out(2) = lwtemp2
            info%flag = zb01_warn_lw

          ELSE
            ! Delete files
            call delete_files(units, number_files, info, filename)
            if(info%flag.lt.0) return

            info%flag = zb01_err_memory_alloc
            RETURN
          END IF
        END IF

        !print *, "restore from file1"
        call read_from_file(w, size_in(1), dest_copy, &
          units, number_files, file_size_copy, info, filename)

        RETURN

      END IF

      ! --------------------------------------
      ! Copy entries
      ! --------------------------------------
      w(dest_copy(1,1):dest_copy(2,1), dest_copy(1,2):dest_copy(2,2)) =&
        wtemp(1:lwtemp1, 1:lwtemp2)

      ! --------------------------------------
      ! Deallocate wtemp
      ! --------------------------------------
      DEALLOCATE (wtemp,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      RETURN

    ELSE
      ! --------------------------------------
      ! Allocation of temporary array not successful or mode = 1
      ! --------------------------------------
      call write_to_file(w, size_in(1), src_copy, &
        units, number_files, file_size_copy, info, filename)
      if(info%flag.lt.0) return

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
      IF (info%stat>0) THEN
        ! --------------------------------------
        ! Allocation of w to desired size failed
        ! Try size=lwtemp
        ! --------------------------------------

        ALLOCATE (w(lwtemp1,lwtemp2),STAT=info%stat)
        IF (info%stat==0) THEN
          ! --------------------------------------
          ! Reassign lw
          ! --------------------------------------
          size_out(1) = lwtemp1
          size_out(2) = lwtemp2
          info%flag = zb01_warn_lw

        ELSE
          ! --------------------------------------
          ! Allocation of w fialed - Delete files
          ! --------------------------------------
          call delete_files(units, number_files, info, filename)
          if(info%flag.lt.0) return

          info%flag = zb01_err_memory_alloc
          RETURN
        END IF
      END IF

      call read_from_file(w, size_out(1), dest_copy, &
        units, number_files, file_size_copy, info, filename)


    END IF

  END IF

END SUBROUTINE zb01_resize2_double

subroutine write_to_file(w, ldw, src, units, number_files, &
      file_size_copy, info, filename)
  integer(long), intent(in) :: ldw ! leading edge of array w
  real(myreal), dimension(ldw, *), intent(in) :: w
  integer(long), dimension(2,2), intent(in) :: src
  integer, dimension(:), allocatable, intent(out) :: units
  integer, intent(out) :: number_files
  integer(long), intent(in) :: file_size_copy
  type(ZB01_info), intent(inout) :: info
  character(len=*), optional, intent(in) :: filename

  integer :: iolen
  CHARACTER (402) :: filename_no
  integer(long) :: ncol
  integer(long) :: nrow
  integer(long) :: total_size
  integer :: i, u
  CHARACTER (10) :: ci        ! Filename extension
  logical :: ex
  integer(long) :: idx1       ! first index into w
  integer(long) :: idx2       ! second index into w
  integer(long) :: written    ! number of bytes written to current file
  integer(long) :: len        ! length to write to file

  !print *, "write to file ", src(:,1), ",", src(:,2)

  ! --------------------------------------
  ! Work out number of files required
  ! --------------------------------------
  INQUIRE (iolength=iolen) w(1,1)
  nrow = src(2,1) - src(1,1)
  ncol = src(2,2) - src(1,2)
  total_size = nrow*ncol*iolen

  number_files = (total_size-1)/file_size_copy + 1
  info%files_used = number_files

  ! --------------------------------------
  ! Check files do not exist
  ! --------------------------------------
  IF (present(filename)) THEN
    DO i = 1, number_files
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      INQUIRE (file=trim(filename_no),exist=ex, iostat=info%iostat)
      IF (info%iostat/=0) THEN
          info%flag = zb01_err_inquire
        return
      ELSE IF (ex) THEN
        info%flag = zb01_err_filename_exists
        return
      END IF
    END DO
  END IF


  ! --------------------------------------
  ! Allocate array for holding unit numbers
  ! --------------------------------------
  i = number_files
  IF (present(filename)) i = 1
  ALLOCATE (units(i),STAT=info%stat)
  IF (info%stat>0) THEN
    info%flag = zb01_err_memory_alloc
    return
  END IF

  ! --------------------------------------
  ! Find unit numbers
  ! --------------------------------------
  call find_units(units, info%iostat)
  if(info%iostat.ne.0) return

  ! --------------------------------------
  ! Open temporary files
  ! --------------------------------------

  IF (present(filename)) THEN
    DO i = 1, number_files
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no),iostat=info%iostat, &
        err=70,status='new',recl=file_size_copy,form='unformatted', &
        action='readwrite')
      CLOSE (unit=u,iostat=info%iostat,err=80,status='keep')
    END DO
  ELSE
    DO i = 1, number_files
      ! write(6,*) 'unit=',units(i)
      u = units(i)
      OPEN (unit=u,iostat=info%iostat,err=70,status='scratch', &
        recl=file_size_copy,form='unformatted',action='readwrite')
    END DO
  END IF

  ! --------------------------------------
  ! Copy entries
  ! --------------------------------------

  idx1 = src(1,1)
  idx2 = src(1,2)
  written = 0
  files: do i = 1, number_files
    if(present(filename)) then
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no), &
        iostat=info%iostat,err=70,status='old',recl=file_size_copy, &
        form='unformatted',action='readwrite', &
        position='rewind')
    else
      u = units(i)
    endif
    do while(file_size_copy - written .gt. 0)
       len = src(2,1) - idx1 + 1 ! amount to reach end of current column
       len = min(len, (file_size_copy-written)/iolen) ! limit to file size
       WRITE (unit=u,iostat=info%iostat,err=90) w(idx1:idx1+len-1,idx2)
       written = written + len*iolen
       idx1 = idx1 + len
       if(idx1.gt.src(2,1)) then
          idx1 = src(1,1)
          idx2 = idx2 + 1
          if(idx2.gt.src(2,2)) then
            if(present(filename)) &
              CLOSE (unit=u,iostat=info%iostat,err=80,status='keep')
            exit files
          endif
       endif
    end do
    if(present(filename)) &
      CLOSE (unit=u,iostat=info%iostat,err=80,status='keep')
  end do files

  return

  70 info%flag = zb01_err_open
  RETURN

  80 info%flag = zb01_err_close
  RETURN

  90 info%flag = zb01_err_write
  RETURN
end subroutine write_to_file

!
! This subroutine fills the array units with available unit numbers
! on which files can be opened
!
subroutine find_units(units, iost)
   integer, dimension(:), intent(out) :: units
   integer, intent(inout) :: iost

   integer :: i ! element of units we need to find
   integer :: u ! current unit to try
   logical :: ex ! .true. if unit exists
   logical :: open ! .true. if unit is not open

   iost = 0 ! initialise in case we do no loop iterations

   u = 8 ! unit to start with
   DO i = 1, size(units)
     DO u = u, huge(0)
       IF (u==100 .OR. u==101 .OR. u==102) CYCLE
       INQUIRE (unit=u,iostat=iost,exist=ex, opened=open)
       if(iost.ne.0) return
       IF (ex .AND. .NOT. open) THEN
         units(i) = u
         EXIT
       END IF
     END DO
     u = units(i) + 1
   END DO
end subroutine find_units

subroutine read_from_file(w, ldw, dest, units, number_files, &
      file_size_copy, info, filename)
  integer(long), intent(in) :: ldw
  real(myreal), dimension(ldw,*), intent(out) :: w
  integer(long), dimension(2,2), intent(in) :: dest
  integer, dimension(:), intent(in) :: units
  integer, intent(in) :: number_files
  integer(long), intent(in) :: file_size_copy
  type(ZB01_info), intent(inout) :: info
  character(len=*), optional, intent(in) :: filename

  integer :: iolen            ! io length of a single element of w
  CHARACTER (402) :: filename_no
  character (10) :: ci
  integer :: i                ! current file
  integer :: u                ! current unit
  integer(long) :: idx1       ! first index into w
  integer(long) :: idx2       ! second index into w
  integer(long) :: done       ! number of bytes read from current file
  integer(long) :: len        ! length to read from file

  !print *, "read from file ", dest(:,1), ",", dest(:,2)

  INQUIRE (iolength=iolen) w(1,1)

  idx1 = dest(1,1)
  idx2 = dest(1,2)
  done = 0
  files: do i = 1, number_files
    if(present(filename)) then
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no), &
        iostat=info%iostat,err=70,status='old',recl=file_size_copy, &
        form='unformatted',action='readwrite', &
        position='rewind')
    else
      u = units(i)
      rewind(u)
    endif
    do while(file_size_copy - done .gt. 0)
       len = dest(2,1) - idx1 + 1 ! amount to reach end of current column
       len = min(len, (file_size_copy-done)/iolen) ! limit to file size
       READ (unit=u,iostat=info%iostat,err=90) w(idx1:idx1+len-1,idx2)
       done = done + len*iolen
       idx1 = idx1 + len
       if(idx1.gt.dest(2,1)) then
          idx1 = dest(1,1)
          idx2 = idx2 + 1
          if(idx2.gt.dest(2,2)) then
            if(present(filename)) &
              CLOSE (unit=u,iostat=info%iostat,err=80,status='delete')
            exit files
          endif
       endif
    end do
    CLOSE (unit=u,iostat=info%iostat,err=80,status='delete')
  end do files

  return

  70 info%flag = zb01_err_open
  RETURN

  80 info%flag = zb01_err_close
  RETURN

  90 info%flag = zb01_err_read
  RETURN
end subroutine read_from_file

subroutine delete_files(units, number_files, info, filename)
  integer, dimension(:), intent(in) :: units
  integer, intent(in) :: number_files
  type(ZB01_info), intent(inout) :: info
  character(len=*), optional, intent(in) :: filename

  character(402) :: filename_no
  character(10) :: ci
  integer :: i, u

  IF (present(filename)) THEN
    do i = 1, number_files
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no), &
        iostat=info%iostat,err=70,status='old', &
        form='unformatted',action='readwrite')
      CLOSE (unit=u,iostat=info%iostat,err=80, &
        status='delete')
    END DO
  ELSE
    DO i = 1, number_files
      u = units(i)
      CLOSE (unit=u,iostat=info%iostat,err=80, &
        status='delete')
    END DO
  END IF

  RETURN

  70 info%flag = zb01_err_open
  RETURN

  80 info%flag = zb01_err_close
  RETURN

end subroutine delete_files

END MODULE hsl_zb01_double
