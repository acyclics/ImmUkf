#include "hungarian.h"


MatrixXd normalize_cost(MatrixXd& cost) {

	int rank = max(cost.rows(), cost.cols());
	MatrixXd normalized_cost = MatrixXd::Zero(rank, rank);
	normalized_cost.topLeftCorner(cost.rows(), cost.cols()) = cost;
	return normalized_cost;
}

vector<int> hungarian_solve(MatrixXd cost) {
	
	cost = normalize_cost(cost);

	const int INF = std::numeric_limits<int>::max();
	const int mrank = cost.rows();
	MatrixXi assignment = MatrixXi::Zero(mrank, mrank);
	int j, k, l, t, q, unmatched;
	double s;

	vector<int> col_mate = vector<int>(mrank, 0);
	vector<int> unchosen_row = vector<int>(mrank, 0);
	vector<int> row_dec = vector<int>(mrank, 0);
	vector<int> slack_row = vector<int>(mrank, 0);

	vector<int> row_mate = vector<int>(mrank, 0);
	vector<int> parent_row = vector<int>(mrank, 0);
	vector<int> col_inc = vector<int>(mrank, 0);
	vector<int> slack = vector<int>(mrank, 0);

	// Begin subtract column minima in order to start with lots of zeroes 12
	for (l = 0; l < mrank; ++l) {
		s = cost(0, l);
		for (k = 1; k < mrank; k++) {
			s = min(s, cost(k, l));
		}
		if (s != 0) {
			for (k = 0; k < mrank; k++) {
				cost(k, l) -= s;
			}
		}
	}
	// End subtract column minima in order to start with lots of zeroes 12

	// Begin initial state 16
	t = 0;
	for (l = 0; l < mrank; l++) {
		row_mate[l] = -1;
		parent_row[l] = -1;
		col_inc[l] = 0;
		slack[l] = INF;
	}
	for (k = 0; k < mrank; k++) {
		s = cost(k, 0);
		for (l = 1; l < mrank; l++) {
			if (cost(k, l) < s) {
				s = cost(k, l);
			}
		}
		row_dec[k] = s;
		for (l = 0; l < mrank; l++) {
			if (s == cost(k, l) && row_mate[l] < 0) {
				col_mate[k] = l;
				row_mate[l] = k;
				goto row_done;
			}
		}
		col_mate[k] = -1;
		unchosen_row[t++] = k;
	row_done:;
	}
	// End initial state 16

	// Begin Hungarian algorithm 18
	if (t == 0) goto done;
	unmatched = t;
	while (true) {
		q = 0;
		while (true) {
			while (q < t) {
				// Begin explore node q of the forest 19
				{
					k = unchosen_row[q];
					s = row_dec[k];
					for (l = 0; l < mrank; l++)
						if (slack[l]) {
							int del;
							del = cost(k, l) - s + col_inc[l];
							if (del < slack[l]) {
								if (del == 0) {
									if (row_mate[l] < 0) goto breakthru;
									slack[l] = 0;
									parent_row[l] = k;
									unchosen_row[t++] = row_mate[l];
								}
								else {
									slack[l] = del;
									slack_row[l] = k;
								}
							}
						}
				}
				// End explore node q of the forest 19
				q++;
			}

			// Begin introduce a new zero into the matrix 21
			s = INF;
			for (l = 0; l < mrank; l++)
				if (slack[l] && slack[l] < s) s = slack[l];
			for (q = 0; q < t; q++) {
				row_dec[unchosen_row[q]] += s;
			}
			for (l = 0; l < mrank; l++)
				if (slack[l]) {
					slack[l] -= s;
					if (slack[l] == 0) {
						// Begin look at a new zero 22
						k = slack_row[l];
						if (row_mate[l] < 0) {
							for (int j = l + 1; j < mrank; j++) {
								if (slack[j] == 0) {
									col_inc[j] += s;
								}
							}
							goto breakthru;
						}
						else {
							parent_row[l] = k;
							unchosen_row[t++] = row_mate[l];
						}
						// End look at a new zero 22
					}
				}
				else {
					col_inc[l] += s;
				}
			// End introduce a new zero into the matrix 21
		}
	breakthru:
		// Begin update the matching 20
		while (true) {
			int j = col_mate[k];
			col_mate[k] = l;
			row_mate[l] = k;
			if (j < 0) break;
			k = parent_row[j];
			l = j;
		}
		// End update the matching 20
		if (--unmatched == 0) goto done;
		// Begin get ready for another stage 17
		t = 0;
		for (l = 0; l < mrank; l++) {
			parent_row[l] = -1;
			slack[l] = INF;
		}
		for (k = 0; k < mrank; k++)
			if (col_mate[k] < 0) {
				unchosen_row[t++] = k;
			}
		// End get ready for another stage 17
	}
done:

	// Begin doublecheck the solution 23
	for (k = 0; k < mrank; k++) {
		for (l = 0; l < mrank; l++) {
			if (cost(k, l) < row_dec[k] - col_inc[l]) assignment(k, l) = -1;
		}
	}
	for (k = 0; k < mrank; k++) {
		l = col_mate[k];
		if (l < 0 || cost(k, l) != row_dec[k] - col_inc[l]) assignment(k, l) = -1;
	}
	// End doublecheck the solution 23
	// End Hungarian algorithm 18

	vector<int> assignments(mrank);

	for (int i = 0; i < mrank; ++i) {
		//assignment(i, col_mate[i]) = 1;
		if (assignment(i, col_mate[i]) != -1)
			assignments[i] = col_mate[i];
		else
			assignments[i] = -1;
	}
	
	return assignments;
}
