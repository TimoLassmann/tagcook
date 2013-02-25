/*
 
 Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of TagDust.
 
 TagDust is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TagDust is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Delve.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#include "nuc_code.h"

void init_nuc_code()
{
	int i;
	for(i = 0;i < 256;i++){
		nuc_code5[i] = 4;
		//	rev_nuc_code16[i] = 15;
	}

	
	nuc_code5[65] = 0;//A Adenine
	nuc_code5[67] = 1;//C	Cytosine
	nuc_code5[71] = 2;//G	Guanine
	nuc_code5[84] = 3;//T	Thymine
	nuc_code5[85] = 3;//U	Uracil
	/*nuc_code5[82] = 5;//R	Purine (A or G)
	 nuc_code5[89] = 6;//Y	Pyrimidine (C, T, or U)
	 nuc_code5[77] = 7;//M	C or A
	 nuc_code5[75] = 8;//K	T, U, or G
	 nuc_code5[87] = 9;//W	T, U, or A
	 nuc_code5[83] = 10;//S	C or G
	 nuc_code5[66] = 11;//B	C, T, U, or G (not A)
	 nuc_code5[68] = 12;//D	A, T, U, or G (not C)
	 nuc_code5[72] = 13;//H	A, T, U, or C (not G)
	 nuc_code5[86] = 14;//V	A, C, or G (not T, not U)
	 nuc_code5[78] = 15;//N	Any base (A, C, G, T, or U)*/
		
	nuc_code5[65+32] = 0;//A Adenine
	nuc_code5[67+32] = 1;//C	Cytosine
	nuc_code5[71+32] = 2;//G	Guanine
	nuc_code5[84+32] = 3;//T	Thymine
	nuc_code5[85+32] = 3;//U	Uracil
	/*nuc_code5[82+32] = 5;//R	Purine (A or G)
	 nuc_code5[89+32] = 6;//Y	Pyrimidine (C, T, or U)
	 nuc_code5[77+32] = 7;//M	C or A
	 nuc_code5[75+32] = 8;//K	T, U, or G
	 nuc_code5[87+32] = 9;//W	T, U, or A
	 nuc_code5[83+32] = 10;//S	C or G
	 nuc_code5[66+32] = 11;//B	C, T, U, or G (not A)
	 nuc_code5[68+32] = 12;//D	A, T, U, or G (not C)
	 nuc_code5[72+32] = 13;//H	A, T, U, or C (not G)
	 nuc_code5[86+32] = 14;//V	A, C, or G (not T, not U)
	 nuc_code5[78+32] = 15;//N	Any base (A, C, G, T, or U)*/
	
	
}

