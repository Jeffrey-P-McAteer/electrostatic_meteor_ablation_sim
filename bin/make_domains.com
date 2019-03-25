#!/bin/csh
# A script to create $1 domain subdirectories for eppic
set verbose = 1
set i=0
while ($i < $1 )
      set name = domain00$i
      if ($i >= 10) set name = domain0$i
      if ($i >= 100) set name = domain$i
      echo $name
      mkdir $name
      @ i +=1
end

end
