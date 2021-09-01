from src.fp import fp_mul, fp_sqr, fp_exp
def novem_1792_3(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_101      = fp_mul(_10,_11)
	_111      = fp_mul(_10,_101)
	_1001     = fp_mul(_10,_111)
	_1011     = fp_mul(_10,_1001)
	_1101     = fp_mul(_10,_1011)
	_1111     = fp_mul(_10,_1101)
	_10001    = fp_mul(_10,_1111)
	_10011    = fp_mul(_10,_10001)
	_10101    = fp_mul(_10,_10011)
	_10111    = fp_mul(_10,_10101)
	_11001    = fp_mul(_10,_10111)
	_11011    = fp_mul(_10,_11001)
	_11101    = fp_mul(_10,_11011)
	_11111    = fp_mul(_10,_11101)
	_111110   = fp_sqr(_11111)
	_111111   = fp_mul(_1,_111110)
	_1000100  = fp_mul(_101,_111111)
	_11111100 = fp_exp(_111111,4)
	_11111111 = fp_mul(_11,_11111100)
	i23       = fp_mul(_10001,_11111111)
	i24       = fp_sqr(_11111111)
	i26       = fp_mul(fp_sqr(i23),_10101)
	i27       = fp_sqr(i24)
	i28       = fp_sqr(i26)
	i29       = fp_sqr(i27)
	i30       = fp_sqr(i28)
	i31       = fp_sqr(i29)
	i32       = fp_sqr(i30)
	i33       = fp_sqr(i31)
	i34       = fp_sqr(i32)
	i35       = fp_sqr(i33)
	i37       = fp_mul(fp_sqr(i34),_1011)
	i38       = fp_sqr(i35)
	i39       = fp_sqr(i37)
	x16       = fp_mul(fp_sqr(i38),_11111111)
	i42       = fp_sqr(i39)
	i43       = fp_sqr(x16)
	i44       = fp_sqr(i42)
	i45       = fp_sqr(i43)
	i46       = fp_sqr(i44)
	i47       = fp_sqr(i45)
	i49       = fp_mul(fp_sqr(i46),_111)
	i50       = fp_sqr(i47)
	i51       = fp_sqr(i49)
	i52       = fp_sqr(i50)
	i53       = fp_sqr(i51)
	i54       = fp_sqr(i52)
	i55       = fp_sqr(i53)
	i56       = fp_sqr(i54)
	i57       = fp_sqr(i55)
	i58       = fp_sqr(i56)
	i59       = fp_sqr(i57)
	i60       = fp_sqr(i58)
	i61       = fp_sqr(i59)
	i62       = fp_sqr(i60)
	i63       = fp_sqr(i61)
	i64       = fp_sqr(i62)
	i66       = fp_mul(fp_sqr(i63),_11101)
	i67       = fp_sqr(i64)
	i68       = fp_sqr(i66)
	i69       = fp_sqr(i67)
	i70       = fp_sqr(i68)
	i71       = fp_sqr(i69)
	i72       = fp_sqr(i70)
	i73       = fp_sqr(i71)
	i74       = fp_sqr(i72)
	x32       = fp_mul(fp_sqr(i73),x16)
	i91       = fp_exp(fp_exp(fp_mul(fp_exp(fp_mul(fp_exp(i74,16),_11001),2),_1),1),256)
	i110      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111111,i91),1),256),_1001),1),256),_1111)
	i138      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i110,64),_1101),1),32768),_11011),1),32)
	i151      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11011,i138),1),8),_111),1),128),_11001)
	i167      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i151,32),_11101),1),16),_1111),1),32)
	i184      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1111,i167),1),64),_10101),1),256),_11)
	i211      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i184,32),_11),1),1024),_10101),1),1024)
	i226      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11101,i211),1),32),_101),1),128),_1011)
	i245      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i226,64),_11111),1),16),_101),1),128)
	i259      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10001,i245),1),32),_11001),1),64),_11001)
	i275      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i259,32),_11011),1),32),_1111),1),16)
	i288      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_101,i275),1),64),_1101),1),16),_111)
	i313      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i288,128),_11001),1),128),_11),1),512)
	i327      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1001,i313),1),16),_101),1),128),_10001)
	i347      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i327,32),_10101),1),64),_11011),1),128)
	i362      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1111,i347),1),1024),_10111),1),4),_1)
	i381      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i362,128),_111111),1),128),_10101),1),8)
	i400      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111,i381),1),256),_11011),1),256),_1011)
	i419      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i400,32),_1001),1),64),_10001),1),64)
	i437      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10111,i419),1),512),_111111),1),64),_1111)
	i457      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i437,32),_1),1),512),_11011),1),16)
	i471      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_1001,i457),1),128),_11001),1),16),_1011)
	i493      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i471,8),_11),1),1024),_11111111),1),128)
	i511      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_11111,i493),1),256),_11111),1),128),_111111)
	i529      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i511,64),_10111),1),16),_111),1),64)
	i544      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_111,i529),1),32),_101),1),128),_11111)
	i563      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i544,64),_10111),1),8),_101),1),256)
	i582      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10011,i563),1),2048),_10101),1),32),_10001)
	i601      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i582,32),_1111),1),16),_11),1),256)
	i618      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(_10111,i601),1),8),_1),1),2048),_11111111)
	i664      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i618,32),_1111),1),64),_11111),1),8589934592)
	i731      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i664),1),4294967296),x32),1),4294967296),x32)
	i829      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i731,4294967296),x32),1),4294967296),x32),1),4294967296)
	i896      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i829),1),4294967296),x32),1),4294967296),x32)
	i994      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i896,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1061     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i994),1),4294967296),x32),1),4294967296),x32)
	i1159     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1061,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1226     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1159),1),4294967296),x32),1),4294967296),x32)
	i1324     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1226,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1391     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1324),1),4294967296),x32),1),4294967296),x32)
	i1489     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1391,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1556     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1489),1),4294967296),x32),1),4294967296),x32)
	i1654     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1556,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1721     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1654),1),4294967296),x32),1),4294967296),x32)
	i1819     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1721,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1886     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1819),1),4294967296),x32),1),4294967296),x32)
	return       fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1886,4294967296),x32),1),4294967296),x32)
