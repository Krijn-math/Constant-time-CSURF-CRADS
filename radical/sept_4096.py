from src.fp import fp_mul, fp_sqr, fp_exp
def sept_4096(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_110      = fp_sqr(_11)
	_111      = fp_mul(_1,_110)
	_1110     = fp_sqr(_111)
	_11100    = fp_sqr(_1110)
	_100011   = fp_mul(_111,_11100)
	_1000110  = fp_sqr(_100011)
	_1000111  = fp_mul(_1,_1000110)
	_1001101  = fp_mul(_110,_1000111)
	_1001110  = fp_mul(_1,_1001101)
	_1110001  = fp_mul(_100011,_1001110)
	_1111111  = fp_mul(_1110,_1110001)
	_11000101 = fp_mul(_1000110,_1111111)
	i15       = fp_mul(_1001101,_11000101)
	i16       = fp_mul(_11000101,i15)
	i18       = fp_mul(fp_sqr(i16),_1000111)
	i19       = fp_mul(i15,i18)
	i20       = fp_mul(_1001110,i19)
	i21       = fp_mul(i19,i20)
	i22       = fp_mul(i16,i21)
	i23       = fp_mul(i21,i22)
	i24       = fp_mul(i20,i23)
	i25       = fp_mul(i18,i24)
	i26       = fp_mul(_1110001,i25)
	i27       = fp_mul(i25,i26)
	i28       = fp_mul(i23,i27)
	i29       = fp_mul(i24,i28)
	i30       = fp_mul(i22,i29)
	i31       = fp_mul(i28,i30)
	i34       = fp_mul(fp_mul(fp_sqr(i31),i31),i26)
	i37       = fp_mul(fp_mul(fp_sqr(i34),i34),i31)
	i38       = fp_mul(i30,i37)
	i39       = fp_mul(i27,i38)
	i41       = fp_mul(fp_sqr(i39),i29)
	i42       = fp_mul(i37,i41)
	i43       = fp_mul(i34,i42)
	i45       = fp_mul(fp_sqr(i43),i39)
	i46       = fp_mul(i41,i45)
	i47       = fp_mul(i45,i46)
	i48       = fp_mul(i43,i47)
	i49       = fp_mul(i47,i48)
	i50       = fp_mul(i48,i49)
	i51       = fp_mul(i38,i50)
	i52       = fp_mul(i49,i51)
	i53       = fp_mul(i42,i52)
	i55       = fp_mul(fp_sqr(i53),i46)
	i57       = fp_mul(fp_sqr(i55),i52)
	i58       = fp_mul(i50,i57)
	i59       = fp_mul(i53,i58)
	i60       = fp_mul(i51,i59)
	i62       = fp_mul(fp_sqr(i60),i58)
	i63       = fp_mul(i59,i62)
	i64       = fp_mul(i55,i63)
	i65       = fp_mul(i63,i64)
	i66       = fp_mul(i57,i65)
	i69       = fp_mul(fp_exp(i66,4),i62)
	i70       = fp_mul(i66,i69)
	i71       = fp_mul(i65,i70)
	i72       = fp_mul(i70,i71)
	i74       = fp_mul(fp_sqr(i72),i69)
	i75       = fp_mul(i71,i74)
	i77       = fp_mul(fp_sqr(i75),i64)
	i78       = fp_mul(i60,i77)
	i79       = fp_mul(i72,i78)
	i80       = fp_mul(i74,i79)
	i82       = fp_mul(fp_sqr(i80),i77)
	i83       = fp_mul(i75,i82)
	i84       = fp_mul(i78,i83)
	i85       = fp_mul(i80,i84)
	i87       = fp_mul(fp_sqr(i85),i79)
	i88       = fp_mul(i82,i87)
	i89       = fp_mul(i84,i88)
	i90       = fp_mul(i85,i89)
	i91       = fp_mul(i87,i90)
	i93       = fp_mul(fp_sqr(i91),i90)
	i95       = fp_mul(fp_sqr(i93),i83)
	i96       = fp_mul(i88,i95)
	i97       = fp_mul(i93,i96)
	i98       = fp_mul(i95,i97)
	i99       = fp_mul(i97,i98)
	i100      = fp_mul(i91,i99)
	i101      = fp_mul(i89,i100)
	i103      = fp_mul(fp_sqr(i101),i96)
	i104      = fp_mul(i99,i103)
	i105      = fp_mul(i100,i104)
	i106      = fp_mul(i98,i105)
	i107      = fp_mul(i103,i106)
	i109      = fp_mul(fp_sqr(i107),i101)
	i110      = fp_mul(i105,i109)
	i113      = fp_mul(fp_mul(fp_sqr(i110),i110),i106)
	i114      = fp_mul(i110,i113)
	i116      = fp_mul(fp_sqr(i114),i109)
	i117      = fp_mul(i104,i116)
	i118      = fp_mul(i116,i117)
	i119      = fp_mul(i114,i118)
	i120      = fp_mul(i117,i119)
	i121      = fp_mul(i113,i120)
	i122      = fp_mul(i107,i121)
	i124      = fp_mul(fp_sqr(i122),i118)
	i125      = fp_mul(i119,i124)
	i126      = fp_mul(i120,i125)
	i127      = fp_mul(i122,i126)
	i128      = fp_mul(i124,i127)
	i129      = fp_mul(i125,i128)
	i130      = fp_mul(i127,i129)
	i131      = fp_mul(i121,i130)
	i132      = fp_mul(i128,i131)
	i133      = fp_mul(i129,i132)
	i134      = fp_mul(i130,i133)
	i135      = fp_mul(i132,i134)
	i136      = fp_mul(i126,i135)
	i137      = fp_mul(_1111111,i136)
	i138      = fp_mul(i131,i137)
	i139      = fp_mul(i133,i138)
	i140      = fp_mul(i135,i139)
	i141      = fp_mul(i139,i140)
	i142      = fp_mul(i134,i141)
	i144      = fp_mul(fp_sqr(i142),i138)
	i145      = fp_mul(i141,i144)
	i146      = fp_mul(i140,i145)
	i147      = fp_mul(i142,i146)
	i148      = fp_mul(i136,i147)
	i149      = fp_mul(i146,i148)
	i150      = fp_mul(i144,i149)
	i151      = fp_mul(i145,i150)
	i152      = fp_mul(i148,i151)
	i153      = fp_mul(i137,i152)
	i154      = fp_mul(i150,i153)
	i155      = fp_mul(i151,i154)
	i156      = fp_mul(i149,i155)
	i157      = fp_mul(i153,i156)
	i158      = fp_mul(i147,i157)
	i159      = fp_mul(i155,i158)
	i160      = fp_mul(i154,i159)
	i161      = fp_mul(i152,i160)
	i162      = fp_mul(i156,i161)
	i163      = fp_mul(i158,i162)
	i164      = fp_mul(i159,i163)
	i165      = fp_mul(i160,i164)
	i166      = fp_mul(i161,i165)
	i167      = fp_mul(i162,i166)
	i168      = fp_mul(i166,i167)
	i169      = fp_mul(i164,i168)
	i171      = fp_mul(fp_sqr(i169),i167)
	i172      = fp_mul(i163,i171)
	i173      = fp_mul(i157,i172)
	i174      = fp_mul(i165,i173)
	i176      = fp_mul(fp_sqr(i174),i173)
	i177      = fp_mul(i168,i176)
	i178      = fp_mul(i171,i177)
	i179      = fp_mul(i172,i178)
	i180      = fp_mul(i174,i179)
	i181      = fp_mul(i178,i180)
	i183      = fp_mul(fp_sqr(i181),i169)
	i184      = fp_mul(i180,i183)
	i185      = fp_mul(i176,i184)
	i186      = fp_mul(i181,i185)
	i187      = fp_mul(i184,i186)
	i188      = fp_mul(i183,i187)
	i190      = fp_mul(fp_sqr(i188),i179)
	i191      = fp_mul(i187,i190)
	i192      = fp_mul(i188,i191)
	i193      = fp_mul(i177,i192)
	i194      = fp_mul(i185,i193)
	i195      = fp_mul(i186,i194)
	i196      = fp_mul(i190,i195)
	i197      = fp_mul(i191,i196)
	i198      = fp_mul(i196,i197)
	i199      = fp_mul(i197,i198)
	i200      = fp_mul(i193,i199)
	i201      = fp_mul(i194,i200)
	i202      = fp_mul(i198,i201)
	i203      = fp_mul(i192,i202)
	i204      = fp_mul(i195,i203)
	i205      = fp_mul(i201,i204)
	i206      = fp_mul(i199,i205)
	i207      = fp_mul(i200,i206)
	i208      = fp_mul(i202,i207)
	i209      = fp_mul(i205,i208)
	i210      = fp_mul(i203,i209)
	i211      = fp_mul(i204,i210)
	i212      = fp_mul(i206,i211)
	i213      = fp_mul(i207,i212)
	i214      = fp_mul(i208,i213)
	i215      = fp_mul(i209,i214)
	i216      = fp_mul(i213,i215)
	i217      = fp_mul(i210,i216)
	i218      = fp_mul(i216,i217)
	i219      = fp_mul(i212,i218)
	i220      = fp_mul(i217,i219)
	i221      = fp_mul(i215,i220)
	i223      = fp_mul(fp_sqr(i221),i214)
	i224      = fp_mul(i218,i223)
	i225      = fp_mul(i211,i224)
	i226      = fp_mul(i219,i225)
	i227      = fp_mul(i225,i226)
	i229      = fp_mul(fp_sqr(i227),i221)
	i230      = fp_mul(i224,i229)
	i231      = fp_mul(i223,i230)
	i232      = fp_mul(i220,i231)
	i233      = fp_mul(i227,i232)
	i234      = fp_mul(i226,i233)
	i235      = fp_mul(i229,i234)
	i238      = fp_mul(fp_exp(i235,4),i234)
	i239      = fp_mul(i235,i238)
	i240      = fp_mul(i238,i239)
	i241      = fp_mul(i230,i240)
	i242      = fp_mul(i231,i241)
	i243      = fp_mul(i233,i242)
	i244      = fp_mul(i240,i243)
	i245      = fp_mul(i239,i244)
	i246      = fp_mul(i242,i245)
	i247      = fp_mul(i243,i246)
	i248      = fp_mul(i244,i247)
	i249      = fp_mul(i241,i248)
	i250      = fp_mul(i247,i249)
	i251      = fp_mul(i246,i250)
	i252      = fp_mul(i248,i251)
	i253      = fp_mul(i232,i252)
	i254      = fp_mul(i250,i253)
	i256      = fp_mul(fp_sqr(i254),i251)
	i257      = fp_mul(i252,i256)
	i258      = fp_mul(i245,i257)
	i259      = fp_mul(i253,i258)
	i260      = fp_mul(i254,i259)
	i261      = fp_mul(i256,i260)
	i262      = fp_mul(i259,i261)
	i263      = fp_mul(i258,i262)
	i264      = fp_mul(i261,i263)
	i265      = fp_mul(i249,i264)
	i267      = fp_mul(fp_sqr(i265),i260)
	i268      = fp_mul(i262,i267)
	i269      = fp_mul(i263,i268)
	i270      = fp_mul(i267,i269)
	i271      = fp_mul(i264,i270)
	i272      = fp_mul(i268,i271)
	i273      = fp_mul(i270,i272)
	i274      = fp_mul(i272,i273)
	x121      = fp_mul(i257,i274)
	i276      = fp_mul(i271,x121)
	i277      = fp_sqr(i276)
	i280      = fp_mul(fp_mul(fp_sqr(i277),i276),i269)
	i281      = fp_mul(i273,i280)
	i282      = fp_mul(i274,i281)
	i285      = fp_mul(fp_exp(i282,4),i265)
	i286      = fp_mul(i280,i285)
	i287      = fp_mul(i276,i286)
	i288      = fp_sqr(i286)
	i290      = fp_mul(fp_mul(i277,i288),i285)
	i291      = fp_mul(i281,i290)
	x128      = fp_mul(i282,i291)
	i683      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i288,340282366920938463463374607431768211456),i291),1),5444517870735015415413993718908291383296),i287),1),680564733841876926926749214863536422912)
	i942      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i290,i683),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1328     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i942,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i1587     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1328),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1973     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1587,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2232     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1973),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i2618     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2232,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2877     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i2618),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i3263     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2877,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i3522     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i3263),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i3908     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i3522,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i4167     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i3908),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	return       fp_mul(fp_exp(i4167,2658455991569831745807614120560689152),x121)
