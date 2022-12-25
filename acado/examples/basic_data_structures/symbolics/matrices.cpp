/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */



 /**
 *    \file examples/basic_data_structures/symbolics.cpp
 *    \author Boris Houska, Hans Joachim Ferreau, Joris Gillis
 *    \date 2010
 */

#include <time.h>

#include <acado/utils/acado_utils.hpp>
#include <acado/user_interaction/user_interaction.hpp>
#include <acado/symbolic_expression/symbolic_expression.hpp>
#include <acado/function/function.hpp>

using namespace std;

USING_NAMESPACE_ACADO

Expression myProduct(const Expression &x, const Expression &y)
{
	return x * y;
}

int main()
{
	// DEFINE VARIABLES:
	// -----------------------

	DifferentialState x("",2, 2);
	Expression z;

	Function f;

	z.appendRows(x);
	z.appendRows(x);

	f << z;

	cout << "x is " << x.getNumRows() << " by " << x.getNumCols() << endl;
	cout << "z is " << z.getNumRows() << " by " << z.getNumCols() << endl;

	// TEST THE FUNCTION f:
	// --------------------
	EvaluationPoint zz(f);

	DVector xx(4), result;
	xx(0) = 2.0;
	xx(1) = 3.0;
	xx(2) = 4.0;
	xx(3) = 5.0;

	zz.setX(xx);
	result = f.evaluate(zz);

	cout << "Result: " << endl << result << endl;

	return 0;
}



