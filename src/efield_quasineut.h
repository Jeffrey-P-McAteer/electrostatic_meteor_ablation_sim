/* code common to the efield_quasineut line of code */
/* int xshift=nx*subdomain.id_number; */
/* int nxT=nx*nsubdomains; */
/* int nrows=nxT*ny; */
/* if (ndim == 3) nrows *= nz; */

inline int gridToR(int ix,int iy,int iz,int nrows) 
{ return (iz*nx*ny + ix*ny + iy + nrows)%nrows; }

/*
inline int rToGrid(int ir,int nx, int ny)
{
  int idx,ix,iy,iz;
  int grid[NDIM];
  if (ndim == 3) {
    idx = ir - xshift*ny*nz;
    iz = idx/(nx*ny);
    iy = idx%ny;
    ix = (idx/ny)%nx;
    grid = {ix,iy,iz};
  }
  else {
    iz = 0; iy = ir%ny; ix = ir/ny-xshift;
    grid = {ix,iy};
  }

  return grid;
}
*/
