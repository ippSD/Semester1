find . -type f | xargs sed -i.bak 's/\#include "stdafx.h"
#include/\#include "stdafx.h"\n\#include/'
