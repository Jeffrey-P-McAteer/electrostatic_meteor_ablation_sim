; This script compresses flux data by a factor of 4 in all dimensions

$mv ppic3d.i ppic3d_bck.i
@eppic.i
@params_in.pro
;cpu, tpool_nthreads = 1
compress_domain_data,'fluxx0.bin','fluxx0.bin',compress=2
compress_domain_data,'fluxy0.bin','fluxy0.bin',compress=2
compress_domain_data,'fluxz0.bin','fluxz0.bin',compress=2

compress_domain_data,'fluxx1.bin','fluxx1.bin',compress=2
compress_domain_data,'fluxy1.bin','fluxy1.bin',compress=2
compress_domain_data,'fluxz1.bin','fluxz1.bin',compress=2

$sed -e 's/nout_avg = 2/nout_avg = 4/' eppic.i > ppic3d.i
;Add some lines with information that allows:analyyzefb3d to work
$echo ";Lines added to allow analyyzefb3d to work" >> ppic3d.i
$echo "nx=nx*nsubdomains" >> ppic3d.i
$echo "nsubdomains=1" >> ppic3d.i





