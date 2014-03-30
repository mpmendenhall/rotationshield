/* 
 * Color.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef COLOR_HH
#define COLOR_HH 1

#include "Vec.hh"
#include <math.h>

template<typename T>
Vec<3,T> rgb2hsv(Vec<3,T> rgb) {
	Vec<3,T> hsv = Vec<3,T>();
	double v = rgb[0];
	double d = rgb[0];
	if(rgb[1]>v) v=rgb[1];
	if(rgb[2]>v) v=rgb[2];
	if(v==0) return hsv;
	if(rgb[1]<d) d=rgb[1];
	if(rgb[2]<d) d=rgb[2];
	d = v-d;
	if(d==0) { hsv[2] = v; return hsv; }
	hsv[1] = d/v;
	if(v == rgb[0])
		hsv[0] = (rgb[1] - rgb[2]) / d;
	else if(v == rgb[1])
		hsv[0] = 2 + (rgb[2] - rgb[0]) / d;
	else
		hsv[0] = 4 + (rgb[0] - rgb[1]) / d;
	hsv[0] *= M_PI/3.0;
	if(hsv[0] < 0)
		hsv[0] += 2*M_PI;
	
	return hsv;
}

template<typename T>
Vec<3,T> hsv2rgb(Vec<3,T> hsv) {
	Vec<3,T> rgb = Vec<3,T>();
	if(hsv[1]==0)
		rgb[0]=rgb[1]=rgb[2]=hsv[2];
	else
	{
		if(hsv[0] < 0) hsv[0] += 2*M_PI;
		double var_h = hsv[0]*3.0/M_PI;
		if(var_h == 6) var_h = 0;
		int var_i = floor(var_h);
		double var_1 = hsv[2] * ( 1 - hsv[1] );
		double var_2 = hsv[2] * ( 1 - hsv[1] * ( var_h - var_i ) );
		double var_3 = hsv[2] * ( 1 - hsv[1] * ( 1 - ( var_h - var_i ) ) );
		
		if(var_i == 0) {rgb[0] = hsv[2]; rgb[1] = var_3; rgb[2] = var_1; }
		else if(var_i == 1) {rgb[0] = var_2; rgb[1] = hsv[2]; rgb[2] = var_1; }
		else if(var_i == 2) {rgb[0] = var_1; rgb[1] = hsv[2]; rgb[2] = var_3; }
		else if(var_i == 3) {rgb[0] = var_1; rgb[1] = var_2; rgb[2] = hsv[2]; }
		else if(var_i == 4) {rgb[0] = var_3; rgb[1] = var_1; rgb[2] = hsv[2]; }
		else {rgb[0] = hsv[2]; rgb[1] = var_1; rgb[2] = var_2; }
	}
	return rgb;
}

#endif
