from src.fp import fp_mul, fp_sqr, fp_exp
def novem_4096(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_110      = fp_sqr(_11)
	_1100     = fp_sqr(_110)
	_1101     = fp_mul(_1,_1100)
	_1110     = fp_mul(_1,_1101)
	_10100    = fp_mul(_110,_1110)
	_10101    = fp_mul(_1,_10100)
	_100011   = fp_mul(_1110,_10101)
	_110001   = fp_mul(_1110,_100011)
	_110100   = fp_mul(_11,_110001)
	_1101000  = fp_sqr(_110100)
	_1111100  = fp_mul(_10100,_1101000)
	_10110000 = fp_mul(_110100,_1111100)
	_10111101 = fp_mul(_1101,_10110000)
	_11100000 = fp_mul(_100011,_10111101)
	i17       = fp_mul(_1111100,_11100000)
	i18       = fp_mul(_110001,i17)
	i20       = fp_mul(fp_sqr(i18),_10110000)
	i21       = fp_mul(_11100000,i20)
	i22       = fp_mul(i18,i21)
	i24       = fp_mul(fp_sqr(i22),i21)
	i27       = fp_mul(fp_mul(fp_sqr(i24),i24)i17)
	i28       = fp_mul(i24,i27)
	i29       = fp_mul(i20,i28)
	i30       = fp_mul(i27,i29)
	i31       = fp_mul(_10111101,i30)
	i33       = fp_mul(fp_sqr(i31),i29)
	i34       = fp_mul(i30,i33)
	i35       = fp_mul(i22,i34)
	i40       = fp_mul(fp_exp(fp_mul(fp_exp(i35,4),i35),2),i31)
	i41       = fp_mul(i28,i40)
	i42       = fp_mul(i35,i41)
	i43       = fp_mul(i33,i42)
	i44       = fp_mul(i40,i43)
	i45       = fp_mul(i34,i44)
	i48       = fp_mul(fp_exp(i45,4),i43)
	i49       = fp_mul(i41,i48)
	i50       = fp_mul(i45,i49)
	i51       = fp_mul(i44,i50)
	i52       = fp_mul(i48,i51)
	i53       = fp_mul(i49,i52)
	i54       = fp_mul(i42,i53)
	i55       = fp_mul(i50,i54)
	i56       = fp_mul(i52,i55)
	i57       = fp_mul(i53,i56)
	i58       = fp_mul(i55,i57)
	i59       = fp_mul(i54,i58)
	i63       = fp_mul(fp_exp(fp_mul(fp_sqr(i59),i59),2),i58)
	i64       = fp_mul(i51,i63)
	i65       = fp_mul(i56,i64)
	i69       = fp_mul(fp_exp(fp_mul(fp_sqr(i65),i65),2),i59)
	i70       = fp_mul(i63,i69)
	i71       = fp_mul(i57,i70)
	i72       = fp_mul(i64,i71)
	i75       = fp_mul(fp_mul(fp_sqr(i72),i72)i65)
	i76       = fp_mul(i69,i75)
	i77       = fp_mul(i71,i76)
	i78       = fp_mul(i70,i77)
	i79       = fp_mul(i72,i78)
	i80       = fp_mul(i76,i79)
	i82       = fp_mul(fp_sqr(i80),i78)
	i84       = fp_mul(fp_sqr(i82),i80)
	i85       = fp_mul(i79,i84)
	i86       = fp_mul(i84,i85)
	i87       = fp_mul(i82,i86)
	i88       = fp_mul(i85,i87)
	i89       = fp_mul(i75,i88)
	i90       = fp_mul(i86,i89)
	i92       = fp_mul(fp_sqr(i90),i88)
	i93       = fp_mul(i87,i92)
	i94       = fp_mul(i92,i93)
	i95       = fp_mul(i90,i94)
	i97       = fp_mul(fp_sqr(i95),i89)
	i98       = fp_mul(i95,i97)
	i99       = fp_mul(i77,i98)
	i100      = fp_mul(i93,i99)
	i101      = fp_mul(i98,i100)
	i102      = fp_mul(i94,i101)
	i103      = fp_mul(i97,i102)
	i106      = fp_mul(fp_mul(fp_sqr(i103),i103)i99)
	i107      = fp_mul(i103,i106)
	i108      = fp_mul(i101,i107)
	i109      = fp_mul(i106,i108)
	i110      = fp_mul(i100,i109)
	i111      = fp_mul(i108,i110)
	i112      = fp_mul(i110,i111)
	i113      = fp_mul(i111,i112)
	i115      = fp_mul(fp_sqr(i113),i109)
	i116      = fp_mul(i107,i115)
	i117      = fp_mul(i112,i116)
	i118      = fp_mul(i113,i117)
	i119      = fp_mul(i102,i118)
	i121      = fp_mul(fp_sqr(i119),i117)
	i122      = fp_mul(i115,i121)
	i123      = fp_mul(i116,i122)
	i124      = fp_mul(i118,i123)
	i125      = fp_mul(i123,i124)
	i126      = fp_mul(i124,i125)
	i127      = fp_mul(i125,i126)
	i128      = fp_mul(i126,i127)
	i129      = fp_mul(i119,i128)
	i130      = fp_mul(i122,i129)
	i131      = fp_mul(i121,i130)
	i133      = fp_mul(fp_sqr(i131),i128)
	i134      = fp_mul(i129,i133)
	i135      = fp_mul(i131,i134)
	i136      = fp_mul(i130,i135)
	i137      = fp_mul(i127,i136)
	i138      = fp_mul(i135,i137)
	i139      = fp_mul(i134,i138)
	i140      = fp_mul(i137,i139)
	i141      = fp_mul(i139,i140)
	i142      = fp_mul(i133,i141)
	i143      = fp_mul(i136,i142)
	i144      = fp_mul(i138,i143)
	i145      = fp_mul(i142,i144)
	i146      = fp_mul(i144,i145)
	i147      = fp_mul(i141,i146)
	i148      = fp_mul(i140,i147)
	i149      = fp_mul(i147,i148)
	i150      = fp_mul(i146,i149)
	i151      = fp_mul(i143,i150)
	i152      = fp_mul(i148,i151)
	i153      = fp_mul(i151,i152)
	i154      = fp_mul(i150,i153)
	i155      = fp_mul(_10101,i154)
	i156      = fp_mul(i153,i155)
	i157      = fp_mul(i154,i156)
	i158      = fp_mul(i145,i157)
	i159      = fp_mul(i152,i158)
	i160      = fp_mul(i157,i159)
	i161      = fp_mul(i158,i160)
	i162      = fp_mul(i156,i161)
	i163      = fp_mul(i159,i162)
	i164      = fp_mul(i149,i163)
	i165      = fp_mul(i162,i164)
	i166      = fp_mul(i155,i165)
	i167      = fp_mul(i163,i166)
	i168      = fp_mul(i164,i167)
	i169      = fp_mul(i166,i168)
	i170      = fp_mul(i161,i169)
	i171      = fp_mul(i160,i170)
	i172      = fp_mul(i170,i171)
	i173      = fp_mul(i171,i172)
	i174      = fp_mul(i168,i173)
	i175      = fp_mul(i173,i174)
	i176      = fp_mul(i169,i175)
	i177      = fp_mul(i167,i176)
	i178      = fp_mul(i165,i177)
	i179      = fp_mul(i174,i178)
	i181      = fp_mul(fp_sqr(i179),i172)
	i182      = fp_mul(i177,i181)
	i184      = fp_mul(fp_sqr(i182),i176)
	i185      = fp_mul(i179,i184)
	i186      = fp_mul(i178,i185)
	i187      = fp_mul(i181,i186)
	i188      = fp_mul(i182,i187)
	i189      = fp_mul(i184,i188)
	i190      = fp_mul(i187,i189)
	i191      = fp_mul(i186,i190)
	i193      = fp_mul(fp_sqr(i191),i175)
	i194      = fp_mul(i185,i193)
	i195      = fp_mul(i193,i194)
	i198      = fp_mul(fp_exp(i195,4),i190)
	i199      = fp_mul(i189,i198)
	i200      = fp_mul(i191,i199)
	i201      = fp_mul(i188,i200)
	i202      = fp_mul(i199,i201)
	i203      = fp_mul(i200,i202)
	i204      = fp_mul(i201,i203)
	i205      = fp_mul(i195,i204)
	i206      = fp_mul(i202,i205)
	i207      = fp_mul(i194,i206)
	i209      = fp_mul(fp_sqr(i207),i203)
	i210      = fp_mul(i205,i209)
	i211      = fp_mul(i198,i210)
	i212      = fp_mul(i207,i211)
	i213      = fp_mul(i211,i212)
	i214      = fp_mul(i206,i213)
	i215      = fp_mul(i212,i214)
	i216      = fp_mul(i214,i215)
	i217      = fp_mul(i209,i216)
	i218      = fp_mul(i216,i217)
	i219      = fp_mul(i204,i218)
	i220      = fp_mul(i218,i219)
	i221      = fp_mul(i219,i220)
	i222      = fp_mul(i220,i221)
	i223      = fp_mul(i221,i222)
	i224      = fp_mul(i217,i223)
	i225      = fp_mul(i210,i224)
	i226      = fp_mul(i215,i225)
	i228      = fp_mul(fp_sqr(i226),i222)
	i229      = fp_mul(i213,i228)
	i230      = fp_mul(i226,i229)
	i231      = fp_mul(i224,i230)
	i232      = fp_mul(i229,i231)
	i233      = fp_mul(i225,i232)
	i234      = fp_mul(i228,i233)
	i235      = fp_mul(i223,i234)
	i236      = fp_mul(i231,i235)
	i237      = fp_mul(i235,i236)
	i238      = fp_mul(i236,i237)
	i239      = fp_mul(i230,i238)
	i240      = fp_mul(i234,i239)
	i241      = fp_mul(i232,i240)
	i242      = fp_mul(i237,i241)
	i243      = fp_mul(i238,i242)
	i245      = fp_mul(fp_sqr(i243),i233)
	i246      = fp_mul(i242,i245)
	i247      = fp_mul(i239,i246)
	i248      = fp_mul(i246,i247)
	i249      = fp_mul(i243,i248)
	i250      = fp_mul(i240,i249)
	i252      = fp_mul(fp_sqr(i250),i247)
	i253      = fp_mul(i241,i252)
	i255      = fp_mul(fp_sqr(i253),i250)
	i256      = fp_mul(i245,i255)
	i257      = fp_mul(i249,i256)
	i259      = fp_mul(fp_sqr(i257),i252)
	i260      = fp_mul(i257,i259)
	i261      = fp_mul(i255,i260)
	i262      = fp_mul(i259,i261)
	i263      = fp_mul(i260,i262)
	i264      = fp_mul(i253,i263)
	i265      = fp_mul(i248,i264)
	i266      = fp_mul(i256,i265)
	i267      = fp_mul(i265,i266)
	i268      = fp_mul(i263,i267)
	x122      = fp_mul(fp_sqr(i268),i268)
	i271      = fp_mul(i261,x122)
	i272      = fp_mul(i262,i271)
	i273      = fp_mul(i264,i272)
	i274      = fp_mul(i267,i273)
	i275      = fp_mul(i273,i274)
	i277      = fp_mul(fp_sqr(i275),i274)
	i278      = fp_mul(i266,i277)
	i279      = fp_mul(i275,i278)
	i281      = fp_mul(fp_sqr(i279),i271)
	i282      = fp_mul(i272,i281)
	i283      = fp_mul(i279,i282)
	i284      = fp_mul(i278,i283)
	i672      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i284,170141183460469231731687303715884105728),i283),1),1361129467683753853853498429727072845824),i281),1),680564733841876926926749214863536422912)
	i931      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_mul(i284,i672),i277),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456)
	i1190     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i282,i931),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282)
	i1576     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1190,340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456)
	i1835     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i282,i1576),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282)
	i2221     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1835,340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456)
	i2480     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i282,i2221),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282)
	i2866     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2480,340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456)
	i3125     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i282,i2866),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282)
	i3511     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i3125,340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456)
	i3770     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i282,i3511),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282)
	i4156     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i3770,340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456),i282),1),340282366920938463463374607431768211456)
	return       fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i282,i4156),1),5316911983139663491615228241121378304),i268),1),4)
