        subroutine savetonetcdf(nrow,nzero,non_row, nzero_id, nzero_value)
          use netcdf
          implicit none
          integer,intent(in) :: nrow,nzero
          integer,intent(in) :: non_row(nrow),nzero_id(nzero)
          real,intent(inout) :: nzero_value(nzero)

          integer :: file_id,non_row_id,nzero_id_id,nzero_value_id,nrow_id,nvalue_id
          integer :: ierr
        
          print*,nrow,nzero
          call err( nf90_create(path='frechetsurf.nc',cmode=NF90_CLOBBER,ncid=file_id))
          call err( nf90_def_dim(file_id,'nrow',nrow,nrow_id))
          call err( nf90_def_dim(file_id,'nzero',nzero,nvalue_id))

          call err( nf90_def_var(file_id,'Non_row',NF90_INT,nrow_id,non_row_id))
          call err( nf90_def_var(file_id,'Non_id',NF90_INT,nvalue_id,nzero_id_id))
          call err( nf90_def_var(file_id,'Non_value',NF90_REAL,nvalue_id,nzero_value_id))

          call err( nf90_enddef(file_id))

          call err( nf90_put_var(file_id,non_row_id,non_row))
          call err( nf90_put_var(file_id,nzero_id_id,nzero_id))
          call err( nf90_put_var(file_id,nzero_value_id,nzero_value))
        
          call err( nf90_close(file_id))
          return
        end subroutine

        subroutine err(jstatus)

        use netcdf

        integer :: jstatus

        if (jstatus .ne. nf90_noerr) then
                print *, trim(nf90_strerror(jstatus))
                stop
        end if

        end subroutine err
