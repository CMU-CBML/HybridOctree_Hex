#ifndef STATICVARS_H
#define STATICVARS_H

static int levelRes[] = {
	1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024
};

static int levelId[] = {
	0,
	1,
	1 + 8,
	1 + 8 + 64,
	1 + 8 + 64 + 512,
	1 + 8 + 64 + 512 + 4096,
	1 + 8 + 64 + 512 + 4096 + 32768,
	1 + 8 + 64 + 512 + 4096 + 32768 + 262144,
	1 + 8 + 64 + 512 + 4096 + 32768 + 262144 + 2097152,
	1 + 8 + 64 + 512 + 4096 + 32768 + 262144 + 2097152 + 16777216,
	1 + 8 + 64 + 512 + 4096 + 32768 + 262144 + 2097152 + 16777216 + 134217728,
	1 + 8 + 64 + 512 + 4096 + 32768 + 262144 + 2097152 + 16777216 + 134217728 + 1073741824
};

static int fId[6][4] = {
	{0, 1, 2, 3},
	{0, 1, 5, 4},
	{0, 3, 7, 4},
	{1, 2, 6, 5},
	{3, 2, 6, 7},
	{4, 5, 6, 7}
};

static int fIdC[6][4] = {
	{0, 1, 2, 3},
	{4, 5, 1, 0},
	{4, 0, 3, 7},
	{5, 6, 2, 1},
	{6, 7, 3, 2},
	{7, 6, 5, 4}
};

static int t1Id[2][13][8] = {
	{{0, 1, 5, 4, 18, 16, 21, 20},
	{1, 2, 6, 5, 16, 17, 22, 21},
	{2, 3, 7, 6, 17, 19, 23, 22},
	{4, 5, 9, 8, 20, 21, 25, 24},
	{5, 6, 10, 9, 21, 22, 26, 25},
	{6, 7, 11, 10, 22, 23, 27, 26},
	{8, 9, 13, 12, 24, 25, 30, 28},
	{9, 10, 14, 13, 25, 26, 31, 30},
	{10, 11, 15, 14, 26, 27, 29, 31},
	{20, 21, 25, 24, 18, 16, 30, 28},
	{22, 23, 27, 26, 17, 19, 29, 31},
	{21, 22, 26, 25, 16, 17, 31, 30},
	{16, 17, 31, 30, 18, 19, 29, 28}},
	{{21, 20, 18, 16, 5, 4, 0, 1},
	{22, 21, 16, 17, 6, 5, 1, 2},
	{23, 22, 17, 19, 7, 6, 2, 3},
	{25, 24, 20, 21, 9, 8, 4, 5},
	{26, 25, 21, 22, 10, 9, 5, 6},
	{27, 26, 22, 23, 11, 10, 6, 7},
	{30, 28, 24, 25, 13, 12, 8, 9},
	{31, 30, 25, 26, 14, 13, 9, 10},
	{29, 31, 26, 27, 15, 14, 10, 11},
	{30, 28, 18, 16, 25, 24, 20, 21},
	{29, 31, 17, 19, 27, 26, 22, 23},
	{31, 30, 16, 17, 26, 25, 21, 22},
	{29, 28, 18, 19, 31, 30, 16, 17}}
};

static int t2Id[2][4][8] = {
	{{0, 8, 9, 2, 1, 12, 13, 3},
	{2, 9, 10, 4, 3, 13, 14, 5},
	{3, 13, 14, 5, 1, 12, 15, 6},
	{4, 10, 11, 7, 5, 14, 15, 6}},
	{{13, 3, 1, 12, 9, 2, 0, 8},
	{14, 5, 3, 13, 10, 4, 2, 9},
	{15, 6, 1, 12, 14, 5, 3, 13},
	{15, 6, 5, 14, 11, 7, 4, 10}}
};

static int t22Id[2][4][8] = {
	{{8, 0, 2, 9, 14, 1, 3, 12},
	{9, 2, 4, 10, 12, 3, 5, 13},
	{12, 3, 5, 13, 14, 1, 6, 15},
	{10, 4, 7, 11, 13, 5, 6, 15}},
	{{3, 12, 14, 1, 2, 9, 8, 0},
	{5, 13, 12, 3, 4, 10, 9, 2},
	{6, 15, 14, 1, 5, 13, 12, 3},
	{6, 15, 13, 5, 7, 11, 10, 4}}
};

static int t3Id[2][3][8] = {
	{{0, 8, 9, 2, 1, 12, 13, 3},
	{2, 9, 10, 4, 3, 13, 14, 5},
	{4, 10, 11, 6, 5, 14, 15, 7}},
	{{13, 3, 1, 12, 9, 2, 0, 8},
	{14, 5, 3, 13, 10, 4, 2, 9},
	{15, 7, 5, 14, 11, 6, 4, 10}}
};

static int t32Id[2][3][8] = {
	{{0, 2, 9, 8, 1, 3, 12, 14},
	{2, 4, 10, 9, 3, 5, 13, 12},
	{4, 6, 11, 10, 5, 7, 15, 13}},
	{{12, 14, 1, 3, 9, 8, 0, 2},
	{13, 12, 3, 5, 10, 9, 2, 4},
	{15, 13, 5, 7, 11, 10, 4, 6}}
};

static int t33Id[2][3][8] = {
	{{8, 0, 2, 9, 12, 1, 3, 13},
	{9, 2, 4, 10, 13, 3, 5, 14},
	{10, 4, 6, 11, 14, 5, 7, 15}},
	{{3, 13, 12, 1, 2, 9, 8, 0},
	{5, 14, 13, 3, 4, 10, 9, 2},
	{7, 15, 14, 5, 6, 11, 10, 4}}
};

static int t34Id[2][3][8] = {
	{{8, 9, 2, 0, 12, 14, 3, 1},
	{9, 10, 4, 2, 14, 15, 5, 3},
	{10, 11, 6, 4, 15, 13, 7, 5}},
	{{3, 1, 12, 14, 2, 0, 8, 9},
	{5, 3, 14, 15, 4, 2, 9, 10},
	{7, 5, 15, 13, 6, 4, 10, 11}}
};

static int t4Id[2][5][8] = {
	{{0, 8, 9, 2, 1, 12, 13, 3},
	{2, 9, 10, 4, 3, 13, 14, 5},
	{4, 10, 11, 6, 5, 14, 15, 7},
	{0, 2, 4, 6, 1, 3, 5, 7},
	{3, 13, 14, 5, 1, 12, 15, 7}},
	{{13, 3, 1, 12, 9, 2, 0, 8},
	{14, 5, 3, 13, 10, 4, 2, 9},
	{15, 7, 5, 14, 11, 6, 4, 10},
	{5, 7, 1, 3, 4, 6, 0, 2},
	{15, 7, 1, 12, 14, 5, 3, 13}}
};

static int t42Id[2][5][8] = {
	{{0, 2, 9, 8, 1, 3, 12, 14},
	{2, 4, 10, 9, 3, 5, 13, 12},
	{4, 6, 11, 10, 5, 7, 15, 13},
	{3, 5, 13, 12, 1, 7, 15, 14},
	{0, 6, 4, 2, 1, 7, 5, 3}},
	{{12, 14, 1, 3, 9, 8, 0, 2},
	{13, 12, 3, 5, 10, 9, 2, 4},
	{15, 13, 5, 7, 11, 10, 4, 6},
	{15, 14, 1, 7, 13, 12, 3, 5},
	{5, 3, 1, 7, 4, 2, 0, 6}}
};

static int t43Id[2][5][8] = {
	{{8, 0, 2, 9, 12, 1, 3, 13},
	{9, 2, 4, 10, 13, 3, 5, 14},
	{10, 4, 6, 11, 14, 5, 7, 15},
	{13, 3, 5, 14, 12, 1, 7, 15},
	{2, 0, 6, 4, 3, 1, 7, 5}},
	{{3, 13, 12, 1, 2, 9, 8, 0},
	{5, 14, 13, 3, 4, 10, 9, 2},
	{7, 15, 14, 5, 6, 11, 10, 4},
	{7, 15, 12, 1, 5, 14, 13, 3},
	{7, 5, 3, 1, 6, 4, 2, 0}}
};

static int t44Id[2][5][8] = {
	{{0, 8, 9, 2, 1, 12, 14, 3},
	{2, 9, 10, 4, 3, 14, 15, 5},
	{4, 10, 11, 6, 5, 15, 13, 7},
	{3, 14, 15, 5, 1, 12, 13, 7},
	{0, 2, 4, 6, 1, 3, 5, 7}},
	{{14, 3, 1, 12, 9, 2, 0, 8},
	{15, 5, 3, 14, 10, 4, 2, 9},
	{13, 7, 5, 15, 11, 6, 4, 10},
	{13, 7, 1, 12, 15, 5, 3, 14},
	{5, 7, 1, 3, 4, 6, 0, 2}}
};

static int pSId[6][2] = {
	{2, 1},
	{2, 0},
	{1, 0},
	{1, 0},
	{2, 0},
	{2, 1}
};

static int pS2Id[6][2] = {
	{3, 4},
	{3, 5},
	{4, 5},
	{4, 5},
	{3, 5},
	{3, 4}
};

static int idTransform[8] = {
	6, 7, 4, 5, 2, 3, 0, 1
};

static int adjP[8][3] = {
	{1, 3, 4},
	{0, 2, 5},
	{1, 3, 6},
	{0, 2, 7},
	{0, 5, 7},
	{1, 4, 6},
	{2, 5, 7},
	{3, 4, 6}
};

static bool towardOutside[8][8][8] = {
	{{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
	{{0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
	{{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
	{{0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0}},
	{{0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0}},
	{{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
	{{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0}},
	{{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}}
};

static double pointOnSurf[146][3] = {
	{ 0.0, 0.0, 1.0 },
	{ 0.3420201433256687, 0.0, 0.9396926207859084 },
	{ 0.3189242855480836, 0.12355192644453342, 0.9396926207859084 },
	{ 0.2527559357867571, 0.23041748059787412, 0.9396926207859084 },
	{ 0.15245149632843955, 0.30616387720913574, 0.9396926207859084 },
	{ 0.03155763752806288, 0.34056114569069434, 0.9396926207859084 },
	{ -0.09359825508738484, 0.32896374433227105, 0.9396926207859084 },
	{ -0.2061131847074455, 0.2729379664507403, 0.9396926207859084 },
	{ -0.29079138662018217, 0.18005040379855497, 0.9396926207859084 },
	{ -0.33619660043916494, 0.06284603641964591, 0.9396926207859084 },
	{ -0.33619660043916494, -0.06284603641964583, 0.9396926207859084 },
	{ -0.29079138662018217, -0.18005040379855503, 0.9396926207859084 },
	{ -0.20611318470744555, -0.2729379664507402, 0.9396926207859084 },
	{ -0.09359825508738492, -0.32896374433227105, 0.9396926207859084 },
	{ 0.03155763752806272, -0.3405611456906944, 0.9396926207859084 },
	{ 0.15245149632843963, -0.30616387720913574, 0.9396926207859084 },
	{ 0.2527559357867571, -0.2304174805978741, 0.9396926207859084 },
	{ 0.3189242855480836, -0.12355192644453346, 0.9396926207859084 },
	{ 0.3420201433256687, 0.0, 0.9396926207859084 },
	{ 0.6427876096865393, 0.0, 0.766044443118978 },
	{ 0.5993815954379041, 0.23220166712762275, 0.766044443118978 },
	{ 0.475025775437305, 0.43304321243580507, 0.766044443118978 },
	{ 0.28651509225520927, 0.5753998723292556, 0.766044443118978 },
	{ 0.05930895822911428, 0.6400455910638803, 0.766044443118978 },
	{ -0.17590717924810528, 0.6182496061102746, 0.766044443118978 },
	{ -0.3873660774325389, 0.512955586012145, 0.766044443118978 },
	{ -0.5465090403901747, 0.3383840716380504, 0.766044443118978 },
	{ -0.6318429291319835, 0.11811191333836743, 0.766044443118978 },
	{ -0.6318429291319835, -0.11811191333836726, 0.766044443118978 },
	{ -0.5465090403901746, -0.33838407163805045, 0.766044443118978 },
	{ -0.38736607743253904, -0.5129555860121447, 0.766044443118978 },
	{ -0.17590717924810542, -0.6182496061102746, 0.766044443118978 },
	{ 0.05930895822911397, -0.6400455910638804, 0.766044443118978 },
	{ 0.28651509225520944, -0.5753998723292556, 0.766044443118978 },
	{ 0.47502577543730506, -0.43304321243580496, 0.766044443118978 },
	{ 0.5993815954379041, -0.23220166712762283, 0.766044443118978 },
	{ 0.6427876096865393, 0.0, 0.766044443118978 },
	{ 0.8660254037844386, 0.0, 0.5000000000000001 },
	{ 0.8075446389876829, 0.31284445982349246, 0.5000000000000001 },
	{ 0.640000495936322, 0.5834375418168269, 0.5000000000000001 },
	{ 0.38602073954358834, 0.7752341508487749, 0.5000000000000001 },
	{ 0.07990674326073385, 0.8623310920878726, 0.5000000000000001 },
	{ -0.2369991014780324, 0.832965440998968, 0.5000000000000001 },
	{ -0.5218969043048337, 0.6911031914823076, 0.5000000000000001 },
	{ -0.736309638274688, 0.4559036264209773, 0.5000000000000001 },
	{ -0.8512796755629922, 0.15913175036229124, 0.5000000000000001 },
	{ -0.8512796755629922, -0.15913175036229102, 0.5000000000000001 },
	{ -0.7363096382746879, -0.4559036264209774, 0.5000000000000001 },
	{ -0.5218969043048339, -0.6911031914823074, 0.5000000000000001 },
	{ -0.23699910147803258, -0.832965440998968, 0.5000000000000001 },
	{ 0.07990674326073345, -0.8623310920878727, 0.5000000000000001 },
	{ 0.3860207395435885, -0.7752341508487748, 0.5000000000000001 },
	{ 0.6400004959363221, -0.5834375418168268, 0.5000000000000001 },
	{ 0.8075446389876829, -0.31284445982349257, 0.5000000000000001 },
	{ 0.8660254037844386, 0.0, 0.5000000000000001 },
	{ 0.984807753012208, 0.0, 0.17364817766693041 },
	{ 0.9183058809859879, 0.3557535935721562, 0.17364817766693041 },
	{ 0.7277817112240621, 0.6634606930336793, 0.17364817766693041 },
	{ 0.43896658858364884, 0.8815637495383915, 0.17364817766693041 },
	{ 0.09086659575717716, 0.9806067367545747, 0.17364817766693041 },
	{ -0.2695054343354901, 0.9472133504425457, 0.17364817766693041 },
	{ -0.5934792621399845, 0.7858935524628853, 0.17364817766693041 },
	{ -0.8373004270103569, 0.5184344754366054, 0.17364817766693041 },
	{ -0.9680395295711485, 0.18095794975801335, 0.17364817766693041 },
	{ -0.9680395295711485, -0.1809579497580131, 0.17364817766693041 },
	{ -0.8373004270103568, -0.5184344754366055, 0.17364817766693041 },
	{ -0.5934792621399846, -0.7858935524628851, 0.17364817766693041 },
	{ -0.26950543433549035, -0.9472133504425457, 0.17364817766693041 },
	{ 0.09086659575717669, -0.9806067367545748, 0.17364817766693041 },
	{ 0.43896658858364906, -0.8815637495383913, 0.17364817766693041 },
	{ 0.7277817112240622, -0.6634606930336792, 0.17364817766693041 },
	{ 0.9183058809859879, -0.3557535935721563, 0.17364817766693041 },
	{ 0.984807753012208, 0.0, 0.17364817766693041 },
	{ 0.984807753012208, 0.0, -0.1736481776669303 },
	{ 0.9183058809859879, 0.3557535935721562, -0.1736481776669303 },
	{ 0.7277817112240621, 0.6634606930336793, -0.1736481776669303 },
	{ 0.43896658858364884, 0.8815637495383915, -0.1736481776669303 },
	{ 0.09086659575717716, 0.9806067367545747, -0.1736481776669303 },
	{ -0.2695054343354901, 0.9472133504425457, -0.1736481776669303 },
	{ -0.5934792621399845, 0.7858935524628853, -0.1736481776669303 },
	{ -0.8373004270103569, 0.5184344754366054, -0.1736481776669303 },
	{ -0.9680395295711485, 0.18095794975801335, -0.1736481776669303 },
	{ -0.9680395295711485, -0.1809579497580131, -0.1736481776669303 },
	{ -0.8373004270103568, -0.5184344754366055, -0.1736481776669303 },
	{ -0.5934792621399846, -0.7858935524628851, -0.1736481776669303 },
	{ -0.26950543433549035, -0.9472133504425457, -0.1736481776669303 },
	{ 0.09086659575717669, -0.9806067367545748, -0.1736481776669303 },
	{ 0.43896658858364906, -0.8815637495383913, -0.1736481776669303 },
	{ 0.7277817112240622, -0.6634606930336792, -0.1736481776669303 },
	{ 0.9183058809859879, -0.3557535935721563, -0.1736481776669303 },
	{ 0.984807753012208, 0.0, -0.1736481776669303 },
	{ 0.8660254037844387, 0.0, -0.4999999999999998 },
	{ 0.807544638987683, 0.3128444598234925, -0.4999999999999998 },
	{ 0.6400004959363221, 0.583437541816827, -0.4999999999999998 },
	{ 0.3860207395435884, 0.775234150848775, -0.4999999999999998 },
	{ 0.07990674326073387, 0.8623310920878727, -0.4999999999999998 },
	{ -0.23699910147803244, 0.8329654409989681, -0.4999999999999998 },
	{ -0.5218969043048338, 0.6911031914823077, -0.4999999999999998 },
	{ -0.7363096382746881, 0.45590362642097737, -0.4999999999999998 },
	{ -0.8512796755629923, 0.15913175036229127, -0.4999999999999998 },
	{ -0.8512796755629923, -0.15913175036229105, -0.4999999999999998 },
	{ -0.736309638274688, -0.4559036264209775, -0.4999999999999998 },
	{ -0.521896904304834, -0.6911031914823075, -0.4999999999999998 },
	{ -0.23699910147803263, -0.8329654409989681, -0.4999999999999998 },
	{ 0.07990674326073345, -0.8623310920878728, -0.4999999999999998 },
	{ 0.38602073954358856, -0.7752341508487749, -0.4999999999999998 },
	{ 0.6400004959363222, -0.5834375418168269, -0.4999999999999998 },
	{ 0.807544638987683, -0.3128444598234926, -0.4999999999999998 },
	{ 0.8660254037844387, 0.0, -0.4999999999999998 },
	{ 0.6427876096865395, 0.0, -0.7660444431189779 },
	{ 0.5993815954379044, 0.23220166712762283, -0.7660444431189779 },
	{ 0.4750257754373052, 0.4330432124358052, -0.7660444431189779 },
	{ 0.2865150922552094, 0.5753998723292558, -0.7660444431189779 },
	{ 0.0593089582291143, 0.6400455910638805, -0.7660444431189779 },
	{ -0.17590717924810534, 0.6182496061102748, -0.7660444431189779 },
	{ -0.38736607743253904, 0.5129555860121451, -0.7660444431189779 },
	{ -0.5465090403901748, 0.3383840716380505, -0.7660444431189779 },
	{ -0.6318429291319837, 0.11811191333836746, -0.7660444431189779 },
	{ -0.6318429291319837, -0.1181119133383673, -0.7660444431189779 },
	{ -0.5465090403901748, -0.33838407163805057, -0.7660444431189779 },
	{ -0.3873660774325392, -0.512955586012145, -0.7660444431189779 },
	{ -0.17590717924810548, -0.6182496061102748, -0.7660444431189779 },
	{ 0.059308958229113994, -0.6400455910638806, -0.7660444431189779 },
	{ 0.2865150922552095, -0.5753998723292557, -0.7660444431189779 },
	{ 0.47502577543730523, -0.4330432124358051, -0.7660444431189779 },
	{ 0.5993815954379044, -0.2322016671276229, -0.7660444431189779 },
	{ 0.6427876096865395, 0.0, -0.7660444431189779 },
	{ 0.3420201433256689, 0.0, -0.9396926207859083 },
	{ 0.3189242855480838, 0.12355192644453347, -0.9396926207859083 },
	{ 0.2527559357867572, 0.23041748059787423, -0.9396926207859083 },
	{ 0.15245149632843963, 0.3061638772091359, -0.9396926207859083 },
	{ 0.03155763752806289, 0.3405611456906945, -0.9396926207859083 },
	{ -0.0935982550873849, 0.3289637443322712, -0.9396926207859083 },
	{ -0.20611318470744558, 0.2729379664507404, -0.9396926207859083 },
	{ -0.29079138662018233, 0.18005040379855508, -0.9396926207859083 },
	{ -0.3361966004391651, 0.06284603641964595, -0.9396926207859083 },
	{ -0.3361966004391651, -0.06284603641964585, -0.9396926207859083 },
	{ -0.2907913866201823, -0.1800504037985551, -0.9396926207859083 },
	{ -0.20611318470744566, -0.27293796645074037, -0.9396926207859083 },
	{ -0.09359825508738497, -0.3289637443322712, -0.9396926207859083 },
	{ 0.03155763752806273, -0.34056114569069457, -0.9396926207859083 },
	{ 0.15245149632843968, -0.30616387720913585, -0.9396926207859083 },
	{ 0.2527559357867572, -0.2304174805978742, -0.9396926207859083 },
	{ 0.3189242855480838, -0.12355192644453351, -0.9396926207859083 },
	{ 0.3420201433256689, 0.0, -0.9396926207859083 },
	{ 0.0, 0.0, -1.0 }
};

static int refinementOnce[8][8] = {
	{19, 1, 4, 0, 5, 6, 9, 8},
	{1, 20, 2, 4, 6, 7, 10, 9},
	{0, 4, 3, 22, 8, 9, 12, 11},
	{4, 2, 21, 3, 9, 10, 13, 12},
	{5, 6, 9, 8, 23, 15, 18, 14},
	{6, 7, 10, 9, 15, 24, 16, 18},
	{8, 9, 12, 11, 14, 18, 17, 26},
	{9, 10, 13, 12, 18, 16, 25, 17}
};
#endif STATICVARS_H