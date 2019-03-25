

#include "timerFile.h"

void reset_timerFile(timerFile &tf) {

  for (timerList::iterator iter=tf.times.begin();
	 iter != tf.times.end();
	 ++iter) {
      iter->second = 0;
    } 
    

}
