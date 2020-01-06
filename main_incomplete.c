#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define GF_size 5
// #define min(a,b) (((a)<(b))? (a) : (b));

int binomialCoeff(int n, int k);
void subset(FILE *fp, int arr[], int data[], int start, int end, int index, int r);
void printsubset(FILE *fp, int arr[], int n, int r);
int printIntersection(FILE *fWriteIntersection, int arr1[], int arr2[], int m, int n);
int printUnion(FILE * fWriteUnion, int arr1[], int arr2[], int m, int n);
int symmDiff(FILE *fWriteDifference, int arr1[], int arr2[], int n, int m);

int main(){
    // Setting
    int N, K, t, d;
	int i, j, k, p, q, ii, r;   // loop index
    // char files[26];
    N = 3;
    K = 4;
	t = 2;
    d = 2;                            // number of copies
    char files[26] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    
    // Prefetching
    int *prefetching_matrix;
    int M, no_present_sym;
    M =  t *((N-1)*t+K-N)/K/(K-1) * binomialCoeff(K, t);
    no_present_sym = N * binomialCoeff(K-1, t-1);
    prefetching_matrix = (int *)malloc(d * no_present_sym * d * no_present_sym * sizeof(int));

    if (prefetching_matrix == 0){
        printf("ERROR: Out of memory\n");
        return 1;
    }
     
    // Assign random number to prefetcing matrix
    // type of mixing files: A+B+C
    for (i = 0; i < d * binomialCoeff(K-1,t-1); i++){
         for (j = 0; j < no_present_sym; j++){
            if (j % (d*binomialCoeff(K-1,t-1)) == i)
                prefetching_matrix[i+j] = rand() % (GF_size - 1) +1;
            else
                prefetching_matrix[i+j]= 0;
         }
    }
    // type of mixing files and segments
    for (k = binomialCoeff(K-1,t-1); k < M; k++){     // big row
        for (i = 0; i < d; i++){                      // small row
            for (j = 0; j < no_present_sym; j++){
                if (j >= d * binomialCoeff(K-1,t-1) && j % d == i)
                    prefetching_matrix[i + k * d + j] = rand() % (GF_size - 1) + 1;
                else
                    prefetching_matrix[i + k * d + j] = 0;
            }
        }
    }
    // this is only true for N=3
    for (i = d * binomialCoeff(K-1,t-1); i < d * (binomialCoeff(K-1,t-1) + 1); i++){
        for (j = d * (no_present_sym - 1); j < d * no_present_sym; j++){
            prefetching_matrix[i + j] = 0;
        }
    }
    for (i = d * (binomialCoeff(K-1,t-1) + 1); i < d * (binomialCoeff(K-1,t-1) + 2); i++){
        for (j = d * (no_present_sym - 2); j < d * (no_present_sym - 1); j++){
            prefetching_matrix[i+j] = 0;
        }
    }
    // posisions for the delivery symbols
    for (i = d * M; i < d * no_present_sym; i++ ){
        for (j = 0; j < no_present_sym; j++){
            prefetching_matrix[i+j] = 0;
        }
    }
    

    /* ==============Delivery================*/
    // Input demands
    char *demand;   
    demand = (char *) malloc(K * sizeof(char));  //memory allocated using malloc
    if(demand == NULL){
        printf("Error! memory not allocated.");
        exit(0);
    }

    printf("Enter demands: ");
    for(i = 0; i < K; ++i)
        scanf("%c", demand + i);
    
    // first, find no. of requests
    int no_requests;
    no_requests = 1;
    for (i = 0; i < K-1; i++){
        if (demand[i+1] != demand[i]){
            no_requests++;
        }
    }   
    printf("The number of requested files is %d\n", no_requests);
    
    // find leaders
    int *leaders;
    leaders = (int *) malloc(no_requests * sizeof(int));
    if(leaders == NULL){
        printf("Error! memory not allocated.");
        exit(0);
    }
    for (j = 0; j < no_requests; j++){
        for(k = 0; k < K; k++){
            if (demand[k] == files[j]){
                leaders[j] = k;
                break;
            }
        }
    }
    printf("The leaders are: ");
    for (j = 0; j< no_requests; j++){
        printf(" %d ", leaders[j]);
    }
    printf("\n");
   
	// find set K(f)
    printf("========== Find Set K(f)=============\n");
    FILE *fWriteKf;
    FILE *fWriteNo;
	int no_requests_this_file;
	char file;
	char filenameKf[sizeof "Koffile*.txt"];
    char filenameNo[sizeof "NoofRequestsofFile*.txt"];
	for (i = 0; i < no_requests; i++){
		file = files[i];
		no_requests_this_file= 0;
		sprintf(filenameKf, "Koffile%c.txt", file);
        sprintf(filenameNo, "NoofRequesteofFile%c.txt", file);
		fWriteKf = fopen(filenameKf, "w");
        fWriteNo = fopen(filenameNo, "w");
		for (j = 0; j < K; j++){
			if (file == demand[j]){
				no_requests_this_file++;
				fprintf(fWriteKf, " %d ",j);
			}
		}
		fprintf(fWriteKf, "\n");
		fprintf(fWriteNo, " %d \n", no_requests_this_file);
		fclose(fWriteKf);
        fclose(fWriteNo);
        printf("Done For File %c \n", file);
	} 

    // set T(i.e. all sets S0) of each file
    printf("============ Find Set S0 =============\n");
    FILE *fWriteS0;
    int T0_size = binomialCoeff(K-1,t);
    char filenameS0[sizeof "S0offile*.txt"];
    int leader;
    int *Kminus1_users;
    for (i = 0; i < no_requests; i++){
        file = files[i];
        sprintf(filenameS0, "S0offile%c.txt", file);
        fWriteS0 = fopen(filenameS0, "w");

        leader = leaders[i];
        Kminus1_users = (int *) malloc((K-1) * sizeof(int));
        if (Kminus1_users == NULL){
            printf("Error! memory not allocated.");
            exit(0);
        }
        j = 0;
        for (k = 0; k < K; k++){
            if (k != leader){
                Kminus1_users[j] = k;
                j++;
            }
        }
        
        printf("Done For File: %c\n", file);
        printsubset(fWriteS0, Kminus1_users, K-1, t); 
        fclose(fWriteS0);
    }

    // Transmission in delivery phase 
    printf("=========== Write Transmission ===========\n");
	FILE *fReadS0;
    FILE *fReadKf;
    FILE *fReadNo;
    FILE *fWriteIntersection;
    FILE *fWriteUnion;
    FILE *fWriteDifference;
    FILE *fWriteDeliveryCopyI;
    FILE *fReadUnion;
    FILE *fReadDifference;
    char filenameIntersection[sizeof "IntersectionSetFile*.txt"];
    char filenameUnion[sizeof "UnionSetFile*.txt"];
    char filenameDifference[sizeof "DifferenceSetFile*.txt"];
    char filenameDeliveryCopyI[sizeof "DeliveryofFile*CopyI.txt"];
    int sizeofintersection;
    int sizeofunion;
    int s[1];
    int sizeofdifference;
    int d_n;   // number of users reqeusting file n
    int uf[1];
    int *S0 = (int *) malloc(t * sizeof(int));
    if (S0 == NULL){
        printf("Error! memory not allocated.");
        exit(0);
    }

    for (i = 0; i < no_requests; i++){   // for each file
        file = files[i];
        uf[0] = leaders[i];
        printf("Doing File %c\n", file);
        sprintf(filenameS0, "S0offile%c.txt", file);
        sprintf(filenameKf, "Koffile%c.txt", file);
        sprintf(filenameNo, "NoofRequesteofFile%c.txt", file);
        fReadS0 = fopen(filenameS0, "r");
        fReadKf = fopen(filenameKf, "r");
        fReadNo = fopen(filenameNo, "r");
        
        // read in each K
        fscanf(fReadNo, "%d", &d_n);   
        int *K = (int *) malloc(d_n * sizeof(int)); // allocate set K(f)
        if (K == NULL){
            printf("Error! memory not allocated.");
            exit(0);
        }
        printf("d_n is %d\n",d_n);
        for (j = 0; j < d_n; j++){    // read in set K(f)
            fscanf(fReadKf, "%d", &K[j]);
            printf("K[%d] is %d\n",j,  K[j]);
        
        } 
        for (k = 0; k < T0_size; k++){     // read in each S0
            printf("Set S0 is:");
            for (j = 0; j < t; j++){
                fscanf(fReadS0, "%d", &S0[j]);
                printf(" %d ", S0[j]);
            }
            printf("\n");

            // S0 \cap K(f)
            sprintf(filenameIntersection, "IntersectionSetFile%c.txt", file);
            fWriteIntersection = fopen(filenameIntersection, "w");
            sizeofintersection = printIntersection(fWriteIntersection, S0 , K, t, d_n);
            
            // uf \cup S0
            sprintf(filenameUnion, "UnionSetFile%c.txt", file);
            fWriteUnion = fopen(filenameUnion, "w");
            sizeofunion = printUnion(fWriteUnion, uf, S0, 1, t);

            // delivery
            sprintf(filenameDeliveryCopyI, "DeliveryofFile%cCopyI.txt", file);
            fWriteDeliveryCopyI = fopen(filenameDeliveryCopyI, "a+");
            fprintf(fWriteDeliveryCopyI, "%c ", file);
            for (j = 0; j < t; j++){
                fprintf(fWriteDeliveryCopyI, "%d ", S0[j]);
            }
            fprintf(fWriteDeliveryCopyI, "[1] ");
            if (sizeofintersection == 0){  // if intersection is empty
                printf("For file A, there is no more");
                fprintf(fWriteDeliveryCopyI, "\n");
            } else {     // if intersection is not empty
                printf("For file A, there is more");
                fprintf(fWriteDeliveryCopyI, " + ");
                // uf \cup S0 \backslash s
                for (p = 0; p < sizeofintersection; p++){  // each s from intersection
                    fscanf(fReadS0, "%d", s);  // retrieve s
                    fReadUnion = fopen(filenameUnion, "r");
					int *ufunionS0;
					ufunionS0 = (int *) malloc(sizeofunion * sizeof(int));
                    for (q = 0; q < sizeofunion; q++){
                        fscanf(fReadUnion, "%d", &ufunionS0[q]); // retrieve uf \cup S0
                    }
                    // perform uf\cup S0 \backslash s, and write it to file
                    fWriteDifference = fopen(filenameDifference, "w");
                    sizeofdifference = symmDiff(fWriteDifference, ufunionS0, s, sizeofunion, 1); 
                    fReadDifference = fopen(filenameDifference, "r");
					int *ufunionS0buts;
					ufunionS0buts = (int *) malloc(sizeofdifference * sizeof(int));
                    for (ii = 0; ii < sizeofdifference; ii++){
                        fscanf(fReadDifference, "%d", &ufunionS0buts[ii]);
                    }
                    // xor the left part
                    fprintf(fWriteDeliveryCopyI, "%c", file);
                    for (r = 0; r < t; r++){  
                        fprintf(fWriteDeliveryCopyI, " %d ", ufunionS0buts[r]);
                    }
                    fprintf(fWriteDeliveryCopyI, " + ");
					free(ufunionS0);
					free(ufunionS0buts);
                    fclose(fReadUnion);
                    fclose(fWriteDifference);
                    fclose(fReadDifference);
                }
                fprintf(fWriteDeliveryCopyI, "\n");
            }

            // Write copy 2:
            fprintf(fWriteDeliveryCopyI, "%c ", file);
            for (j = 0; j < t; j++){
                fprintf(fWriteDeliveryCopyI, "%d ", S0[j]);
            }
            fprintf(fWriteDeliveryCopyI, "[1] + %c", file);
            for (j = 0; j < t; j++){
                fprintf(fWriteDeliveryCopyI, "%d ", S0[j]);
            }
            fprintf(fWriteDeliveryCopyI, "[2] ");
            if (sizeofintersection == 0){  // if intersection is empty
                fprintf(fWriteDeliveryCopyI, "\n");
            } else {     // if intersection is not empty
                fprintf(fWriteDeliveryCopyI, " + ");
                // uf \cup S0 \backslash s
                for (p = 0; p < sizeofintersection; p++){  // each s from intersection
                    fscanf(fReadS0, "%d", s);  // retrieve s
                    fReadUnion = fopen(filenameUnion, "r");
					int *ufunionS0;
					ufunionS0 = (int *)malloc(sizeofunion * sizeof(int));
                    for (q = 0; q < sizeofunion; q++){
                        fscanf(fReadUnion, "%d", &ufunionS0[q]); // retrieve uf \cup S0
                    }
                    // perform uf\cup S0 \backslash s, and write it to file
                    fWriteDifference = fopen(filenameDifference, "w");
                    sizeofdifference = symmDiff(fWriteDifference, ufunionS0, s, sizeofunion, 1); 
                    fReadDifference = fopen(filenameDifference, "r");
					int *ufunionS0buts;
					ufunionS0buts = (int *)malloc(sizeofdifference * sizeof(int));
                    for (ii = 0; ii < sizeofdifference; ii++){
                        fscanf(fReadDifference, "%d", &ufunionS0buts[ii]);
                    }
                    // xor the left part
                    fprintf(fWriteDeliveryCopyI, "%c", file);
                    for (r = 0; r < t; r++){  
                        fprintf(fWriteDeliveryCopyI, " %d ", ufunionS0buts[r]);
                    }
                    fprintf(fWriteDeliveryCopyI, "[1] + %c", file);
                    for (r = 0; r < t; r++){  
                        fprintf(fWriteDeliveryCopyI, " %d ", ufunionS0buts[r]);
                    }
                    fprintf(fWriteDeliveryCopyI, "[2]");
					free(ufunionS0);
					free(ufunionS0buts);
                    fclose(fReadUnion);
                    fclose(fWriteDifference);
                    fclose(fReadDifference);
                }
                fprintf(fWriteDeliveryCopyI, "\n");
            } // end if-else

            fclose(fWriteIntersection);
            fclose(fWriteUnion);
            fclose(fWriteDeliveryCopyI);
        } // end for-loop of S0
        free(K);
        fclose(fReadS0);
        fclose(fReadKf);
        fclose(fReadNo);
       // fclose(ReadDifference);
        //fclose(fWriteUnion);
        printf("Done For File %c\n", file);
    } // end for-loop of file
    free(S0);
    
    
    /* ============ Last Step: Free Memory ==============*/
    free(prefetching_matrix);
    free(demand);
    free(leaders);
    return 0;
} // end main function



// Returns value of Binomial Coefficient C(n, k)
int binomialCoeff(int n, int k) {
    int C[n+1][k+1];
    int i, j;

    // Caculate value of Binomial Coefficient in bottom up manner
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= min(i, k); j++){
            // Base Cases
            if (j == 0 || j == i)
                C[i][j] = 1;

            // Calculate value using previosly stored values
            else
                C[i][j] = C[i-1][j-1] + C[i-1][j];
        }
    }

    return C[n][k];
}

// A utility function to return minimum of two integers
// int min(int a, int b){
//    return (a<b)? a: b;
//}

/*
 * // create d-by-d random diagonal matrix
* void create_random_matrix(int *random_matrix, int d){
    for (int i = 1; i <= M; i++)
        for (int j = 1; j <= no_present_sym; j++){
            if (i == j)
                random_matrix(i,j) = rand() % (GF_size - 1) + 1;
            else
                random_matrix(i,j) = 0;    
        }
}
*/
/*  Function to generate subset  */
void subset(FILE *fp, int arr[], int data[], int start, int end, int index, int r){
    int j, i;
    if (index == r) {
        for (j = 0; j < r; j++)
            fprintf(fp, "%d ", data[j]);
        fprintf(fp, "\n");
        return;
    }
    for (i = start; i <= end && end - i + 1 >= r - index; i++) {
        data[index] = arr[i];
        subset(fp, arr, data, i+1, end, index+1, r);
    }
}
/*  End of subset()  */
 
/*  Function to print the subset  */ 
void printsubset(FILE *fp, int arr[], int n, int r){
	int *data;
	data = (int *)malloc(r * sizeof(int));
    subset(fp, arr, data, 0, n - 1, 0, r);
}

/* Function prints union of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
int printUnion(FILE *fWriteUnion, int arr1[], int arr2[], int m, int n){
    int i = 0, j = 0;
    int sizeofunion = 0;
    while (i < m && j < n){
        if (arr1[i] < arr2[j]){
            fprintf(fWriteUnion, " %d ", arr1[i++]);
            sizeofunion++;
        } else if (arr2[j] < arr1[i]){
            fprintf(fWriteUnion, " %d ", arr2[j++]);
            sizeofunion++;
        } else{
            fprintf(fWriteUnion, " %d ", arr2[j++]);
            i++;
            sizeofunion++;
        }
    }

  /* Print remaining elements of the larger array */
    while(i < m){
        fprintf(fWriteUnion, " %d ", arr1[i++]);
        sizeofunion++;
    }
    while(j < n){
        fprintf(fWriteUnion, " %d ", arr2[j++]);
        sizeofunion++;
    }
    return sizeofunion;
}

// two sorted arrays
/* Function prints Intersection of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
int printIntersection(FILE *fWriteIntersection, int arr1[], int arr2[], int m, int n) {
    int i = 0, j = 0;
    int sizeofIntersection = 0;
    while (i < m && j < n){
        if (arr1[i] < arr2[j])
            i++;
        else if (arr2[j] < arr1[i])
            j++;
        else{ /* if arr1[i] == arr2[j] */
            fprintf(fWriteIntersection, " %d ", arr2[j++]);
            i++;
            sizeofIntersection++;
        }
    }
    return sizeofIntersection;
}


int symmDiff(FILE *fWriteDifference, int arr1[], int arr2[], int n, int m) {
    // Traverse both arrays simultaneously.
    int i = 0, j = 0;
    int sizeofdifference = 0;
    while (i < n && j < m) {
        // Print smaller element and move 
        // ahead in array with smaller element
        if (arr1[i] < arr2[j]){
            fprintf(fWriteDifference, "%d", arr1[i]);
            i++;
            sizeofdifference++;
        }
        else if (arr2[j] < arr1[i]){
            fprintf(fWriteDifference, "%d", arr2[j]);
            j++;
            sizeofdifference++;
        }
        // If both elements same, move ahead
        // in both arrays.
        else{
            i++;
            j++;
        }
    }
    return sizeofdifference;
}
