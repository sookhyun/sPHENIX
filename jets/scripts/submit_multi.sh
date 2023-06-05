#!/bin/csh 

#         << " 4. Jet Radius" << endl
#         << " 7. ISR" << endl
#         << " 8. Multiparton Interaction" << endl
#         << " 9. Hadronization" << endl
#         << " 10. sqrtS " << endl
#         << " 11. Omega_q selected " << endl
#         << " 12. File index" << endl
#

set nruns = 1000

set index = 0
while ($index < $nruns)
  ./submit_select.sh $index
  #echo 'submitting job ' $index
  @ index++
end

