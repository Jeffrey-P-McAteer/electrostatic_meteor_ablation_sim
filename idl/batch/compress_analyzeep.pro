; This script compresses a data set by a factor of 4 in all dimensions
; and then analyzes it.


@eppic.i
@params_in.pro
;cpu, tpool_nthreads = 1

file_mkdir, 'compress'
file_mkdir, 'compress/domain000'
compress_domain_data,'phi.bin','compress/domain000/phi.bin',compress=2
compress_domain_data,'den0.bin','compress/domain000/den0.bin',compress=2
compress_domain_data,'den1.bin','compress/domain000/den1.bin',compress=2
file_copy, 'eppic.i','compress/eppic.i',/overwrite
$cp -pf domain000/*.out compress/domain000/

;Modify and add some lines with information that allows analyyzefbep
;to work on a single domain data set
$sed -e 's/nout_avg = 2/nout_avg = 4/' eppic.i > compress/eppic.i
cd,'compress'
$echo ";Lines added to allow analyyzefb3d to work" >> eppic.i
$echo "nx=nx*nsubdomains" >> eppic.i
$echo "nsubdomains=1" >> eppic.i

@analyzefbep



