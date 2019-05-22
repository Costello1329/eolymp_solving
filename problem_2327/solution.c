#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>




void _mergesort (int from, int to, int *in, int *out) {
	if (from == to) {
		return;
	}

	const int mid = (from + to) / 2;
	_mergesort(from, mid, in, out);
	_mergesort(mid + 1, to, in, out);

	int i = from;
	int j = mid + 1;

	for (int k = 0; k < to - from + 1; k ++) {
		if (i <= mid && j <= to) {
			if (in[i] < in[j]) {
				out[k] = in[i];
				i ++;
			}

			else {
				out[k] = in[j];
				j ++;
			}
		}

		else if (j > to) {
			out[k] = in[i];
			i ++;
		}

		else {
			out[k] = in[j];
			j ++;
		}
	}

	for (int i = 0; i < to - from + 1; i ++) {
		in[i + from] = out[i];
	}
}



int main () {
	int n;
	scanf("%d", &n);

	int *in = (int *) malloc(sizeof(int) * n);
	int *out = (int *) malloc(sizeof(int) * n);

	for (int i = 0; i < n; i ++)
		scanf("%d", &in[i]);

	_mergesort(0, n - 1, in, out);

	for (int i = 0; i < n; i ++)
		printf("%d ", in[i]);

	printf("\n");

	free(in);
	free(out);
}

