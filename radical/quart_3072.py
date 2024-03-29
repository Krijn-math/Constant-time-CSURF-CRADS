from src.fp import fp_mul, fp_sqr, fp_exp
def quart_3072(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_100      = fp_mul(_1,_11)
	_101      = fp_mul(_1,_100)
	_111      = fp_mul(_10,_101)
	_1001     = fp_mul(_10,_111)
	_1011     = fp_mul(_10,_1001)
	_1101     = fp_mul(_10,_1011)
	_1111     = fp_mul(_10,_1101)
	_10011    = fp_mul(_100,_1111)
	_10101    = fp_mul(_10,_10011)
	_10111    = fp_mul(_10,_10101)
	_11001    = fp_mul(_10,_10111)
	_11011    = fp_mul(_10,_11001)
	_11101    = fp_mul(_10,_11011)
	_11111    = fp_mul(_10,_11101)
	_100001   = fp_mul(_10,_11111)
	_100101   = fp_mul(_100,_100001)
	_100111   = fp_mul(_10,_100101)
	_101001   = fp_mul(_10,_100111)
	_101011   = fp_mul(_10,_101001)
	_101101   = fp_mul(_10,_101011)
	_101111   = fp_mul(_10,_101101)
	_110001   = fp_mul(_10,_101111)
	_110011   = fp_mul(_10,_110001)
	_110101   = fp_mul(_10,_110011)
	_111001   = fp_mul(_100,_110101)
	_111101   = fp_mul(_100,_111001)
	_111111   = fp_mul(_10,_111101)
	_1111110  = fp_sqr(_111111)
	_1111111  = fp_mul(_1,_1111110)
	_10110100 = fp_mul(_110101,_1111111)
	i53       = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(_10110100,64),_101111),1),64),_1111),1),128)
	i73       = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11101,i53),1),512),_101001),1),256),_10101)
	i96       = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i73,128),_101111),1),512),_100001),1),32)
	i113      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10011,i96),1),64),_1001),1),256),_110001)
	i133      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i113,32),_10101),1),16),_101),1),512)
	i147      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111001,i133),1),16),_101),1),128),_1011)
	i172      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i147,256),_110101),1),32),_10111),1),1024)
	i190      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_100101,i172),1),64),_101011),1),512),_100101)
	i214      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i190,128),_101111),1),64),_11),1),512)
	i234      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11101,i214),1),128),_1011),1),1024),_10011)
	i261      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i234,128),_111),1),1024),_100111),1),256)
	i282      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_100111,i261),1),512),_110101),1),512),_1111111)
	i303      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i282,32),_11111),1),512),_100101),1),32)
	i319      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10111,i303),1),64),_11011),1),128),_111101)
	i340      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i319,64),_110101),1),64),_10111),1),128)
	i361      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111111,i340),1),128),_10011),1),2048),_1001)
	i386      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i361,1024),_111111),1),128),_100111),1),64)
	i400      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11011,i386),1),128),_11011),1),16),_111)
	i423      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i400,1024),_101011),1),64),_101011),1),32)
	i441      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10101,i423),1),64),_11111),1),512),_101001)
	i463      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i441,64),_100111),1),8),_101),1),2048)
	i477      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111101,i463),1),16),_1101),1),128),_1001)
	i502      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i477,512),_11101),1),256),_100101),1),64)
	i518      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_101101,i502),1),64),_100001),1),128),_110001)
	i544      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i518,64),_101),1),128),_101),1),2048)
	i557      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_mul(_100001,i544),2),_1),1),512),_1111111)
	i580      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i557,64),_11001),1),256),_101011),1),128)
	i595      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11101,i580),1),32),_1001),1),128),_11111)
	return       fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i595,128),_110011),1),32),_11101),1),13540631571091461853163616343597829035214122634757923946150980204596801748539633714799749980264179980479764347520065911539820003456901359936565697225764325861376456757098784868801942261231107873715432358778281570209476232615685170429490116466837180248143047642366386324907479955857757773156061115886947408860370347929954110401396697023036971669330349689165406178637931741075554388570447364695854236037800204543434404957923178194336489464779659687483143656222678723617764586826350928747474366746593079282080501565859402166884818402721909226979370809408774598305372366324022526391307845011736946602338267951389450620059880839845574200135281351584330463507370495974012473127571340647450941142988034283084116727091178291918317863076646033235212001226644720063655904142163968)
