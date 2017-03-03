        subroutine savetonetcdf(non_row, nzero_id, nzero_value)
          use netcdf
          implicit none
          integer,intent(in),dimension(:)::non_row,nzero_id
          integer,parameter::dp = SELECTED_REAL_KIND(15,307)
          real(kind=dp),intent(in),dimension(:)::nzero_value

          integer :: file_id,non_row_id,nzero_id_id,nzero_value_id,nrow_id,nvalue_id
          integer :: nrow,nzero
          integer :: ierr
          nrow = size(non_row)
          nzero = size(nzero_id)
          ierr = nf90_create(path='frechet.nc',cmode=NF90_CLOBBER,ncid=file_id)
          ierr = nf90_def_dim(file_id,'nrow',nrow,nrow_id)
          ierr = nf90_def_dim(file_id,'nzero',nzero,nvalue_id)

          ierr = nf90_def_var(file_id,'Non_row',NF90_INT,nrow_id,non_row_id)
          ierr = nf90_def_var(file_id,'Non_id',NF90_INT,nvalue_id,nzero_id_id)
          ierr = nf90_def_var(file_id,'Non_value',NF90_REAL,nvalue_id,nzero_value_id)

          ierr = nf90_enddef(file_id)

          ierr = nf90_put_var(file_id,non_row_id,non_row)
          ierr = nf90_put_var(file_id,nzero_id_id,nzero_id)
          ierr = nf90_put_var(file_id,nzero_value_id,nzero_value)
        
          ierr = nf90_close(file_id)
          return
        end subroutine
