from src.fp import fp_mul, fp_sqr, fp_exp
def sq_2048(x):
	_1 = x
	_10     = fp_sqr(_1)
	_11     = fp_mul(_1,_10)
	_101    = fp_mul(_10,_11)
	_111    = fp_mul(_10,_101)
	_1001   = fp_mul(_10,_111)
	_1011   = fp_mul(_10,_1001)
	_1101   = fp_mul(_10,_1011)
	_1111   = fp_mul(_10,_1101)
	_10001  = fp_mul(_10,_1111)
	_10011  = fp_mul(_10,_10001)
	_10101  = fp_mul(_10,_10011)
	_10111  = fp_mul(_10,_10101)
	_11001  = fp_mul(_10,_10111)
	_11011  = fp_mul(_10,_11001)
	_11101  = fp_mul(_10,_11011)
	_11111  = fp_mul(_10,_11101)
	_110010 = fp_mul(_10011,_11111)
	_111111 = fp_mul(_1101,_110010)
	i36     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(_110010,8),_101),1),128),_11101),1),64)
	i52     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1001,i36),1),128),_111111),1),64),_10011)
	i75     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i52,64),_1011),1),1024),_11111),1),32)
	i91     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11,i75),1),256),_10011),1),32),_11001)
	i110    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i91,4),_11),1),2048),_10111),1),16)
	i123    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1101,i110),1),64),_10001),1),16),_11)
	i143    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i123,128),_1001),1),128),_10011),1),16)
	i158    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111,i143),1),64),_1001),1),64),_11111)
	i177    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i158,32),_1101),1),128),_1011),1),32)
	i192    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111,i177),1),8),_1),1),512),_11111)
	i216    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i192,128),_1011),1),256),_10101),1),128)
	i228    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11101,i216),1),16),_1101),1),32),_1011)
	i249    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i228,64),_11101),1),256),_1111),1),32)
	i271    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_101,i249),1),1024),_10001),1),512),_111111)
	i295    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i271,1024),_111),1),32),_101),1),128)
	i315    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11,i295),1),1024),_1111),1),128),_10001)
	i336    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i315,32),_11011),1),256),_111),1),64)
	i349    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111,i336),1),128),_11001),1),8),_101)
	i369    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i349,32),_11),1),512),_10011),1),16)
	i384    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1011,i369),1),32),_101),1),128),_11111)
	i402    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i384,64),_11101),1),16),_1101),1),64)
	i419    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11111,i402),1),64),_1001),1),256),_10001)
	i442    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i419,8),_11),1),1024),_10011),1),256)
	i454    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10111,i442),1),64),_1111),1),8),_1)
	i471    = fp_exp(fp_exp(fp_mul(fp_exp(fp_mul(fp_exp(i454,512),_111111),2),_1),1),32)
	i488    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1101,i471),1),512),_11111),1),32),_1001)
	i507    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i488,8),_11),1),256),_1001),1),64)
	i519    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11111,i507),1),64),_10011),1),8),_101)
	i540    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i519,128),_11001),1),64),_10001),1),64)
	i553    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10001,i540),1),32),_10101),1),32),_1111)
	i572    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i553,64),_1011),1),16),_111),1),128)
	i586    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10111,i572),1),16),_1011),1),128),_11011)
	i605    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i586,4),_11),1),256),_11101),1),128)
	return     fp_exp(fp_exp(fp_mul(_1001,i605),1),602578106730258147145029151507078528228023311986422237839918759766314347897950800167271256026836512415862191570222001196965604244698739581242403246972693662863803334845406306192695519479420374919211385784346955014399336464476149777365346780524876991249453910335575169407368338820404357102948540991946467592592242277305397985763558002890689612520144896962725248428723070869161657795378765951218343795565126561748209474676492753131460831100308246714251190084829184)
