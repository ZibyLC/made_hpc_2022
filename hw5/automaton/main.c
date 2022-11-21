#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define SIZE 100
#define MAX_ITERS 20

#define DEBUG 1
#define VERBOSE 1
#define CYCLIC_WORLD 1
#define GHOST_OFFSET 1

#define RULE 110

#define FIX_UNUSED(X) (void)(X)
#define CHECK_RETURN(RET, MSG)                                                    \
  if (RET) {                                                                   \
    printf("Error: %s: ret = %d", MSG, RET);                                   \
    return -1;                                                                 \
  }

int get_N(int p_rank, int p_size) {
    FIX_UNUSED(p_rank);
    assert(SIZE > p_size);
    return SIZE / p_size +
           2 * GHOST_OFFSET +
           (p_rank < SIZE % p_size);
}

int prev_rank(int p_rank, int p_size) {
    int prev = p_rank - 1;
    if (CYCLIC_WORLD)
        return (prev + p_size) % p_size;
    else
        return prev;
}

int next_rank(int p_rank, int p_size) {
    int next = p_rank + 1;
    if (CYCLIC_WORLD)
        return next % p_size;
    else
        return next;
}

int begin(int N) {
    FIX_UNUSED(N);
    return GHOST_OFFSET;
}

int end(int N) {
    assert(N >= 3 * GHOST_OFFSET);
    return N - GHOST_OFFSET;
}

int count(int N) {
    return end(N) - begin(N);
}

void init_state(int p_rank, int p_size, unsigned int *state) {
    int N = get_N(p_rank, p_size);
    unsigned int seed = (p_rank + 1) * (unsigned int) time(NULL);
    for (int i = 0; i < begin(N); i++)
        state[i] = 0;
    for (int i = begin(N); i < end(N); i++)
        state[i] = (unsigned int) rand_r(&seed) % 2;
    for (int i = end(N); i < N; i++)
        state[i] = 0;
}

int sent_recv_ghost_cells(int p_rank, int p_size, int to_next, unsigned int *state) {
    int N = get_N(p_rank, p_size);
    int next = next_rank(p_rank, p_size);
    int prev = prev_rank(p_rank, p_size);
    int dst_prank = to_next ? next : prev;
    int src_prank = to_next ? prev : next;
    int src_offset = to_next ? end(N) - GHOST_OFFSET : begin(N);
    int dst_offset = to_next ? 0 : end(N);
    if (dst_prank >= 0 && dst_prank < p_size)
        CHECK_RETURN(
                MPI_Send(state + src_offset, GHOST_OFFSET, MPI_UNSIGNED, dst_prank, 0, MPI_COMM_WORLD),
                "MPI_Send");
    if (src_prank >= 0 && src_prank < p_size)
        CHECK_RETURN(
                MPI_Recv(state + dst_offset, GHOST_OFFSET, MPI_UNSIGNED, src_prank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE),
                "MPI_Recv");
    return 0;
}

int make_step(int ghost_lo, int ghost_hi, int lo, int hi, unsigned int *state) {
    assert(RULE >= 0 && RULE < 256);
    int converged = 1;
    unsigned int prev = ghost_lo;
    for (int i = lo; i < hi; i++) {
        unsigned int curr = state[i];
        unsigned int next = i + 1 < hi ? state[i + 1] : ghost_hi;
        unsigned int val = (next) ^ (curr << 1) ^ (prev << 2);
        state[i] = 0;
        for (unsigned int pos = 0; pos < 8; pos++) {
            if (RULE & (1 << pos))
                state[i] |= (val == pos);
        }
        if (DEBUG) {
            if (RULE == 110)
                assert(state[i] == (val == 1 || val == 2 || val == 3 || val == 5 || val == 6));
            else if (RULE == 90)
                assert(state[i] == (val == 1 || val == 3 || val == 4 || val == 6));
            else if (RULE == 30)
                assert(state[i] == (val == 1 || val == 2 || val == 3 || val == 4));
        }
        if (curr != state[i])
            converged = 0;
        prev = curr;
    }
    return converged;
}

int gather_states(int p_rank, int p_size, const unsigned int *state, unsigned int *full_state) {
    int N = get_N(p_rank, p_size);
    if (p_rank == 0) {
        for (int i = begin(N), j = 0; i < end(N); i++, j++) {
            full_state[j] = state[i];
        }
        for (int src = 1, offset = count(N); src < p_size; src++) {
            int src_N = get_N(src, p_size);
            int src_count = count(src_N);
            CHECK_RETURN(
                    MPI_Recv(full_state + offset, src_count, MPI_UNSIGNED, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE),
                    "MPI_Recv");
            offset += src_count;
        }
    } else {
        int dst = 0;
        CHECK_RETURN(
                MPI_Send(state + begin(N), count(N), MPI_UNSIGNED, dst, 0, MPI_COMM_WORLD),
                "MPI_Send");
    }
    return 0;
}

int run_wolfram_automaton(int p_rank, int p_size) {
    int N = get_N(p_rank, p_size);
    unsigned int *state = NULL, *full_state = NULL, *prev_full_state = NULL;
    state = (unsigned int *) malloc(N * sizeof(unsigned int));
    if (p_rank == 0 && DEBUG) {
        full_state = (unsigned int *) malloc(SIZE * sizeof(unsigned int));
        prev_full_state = (unsigned int *) malloc(SIZE * sizeof(unsigned int));
    }

    // initialize cells
    init_state(p_rank, p_size, state);
    if (DEBUG) {
        gather_states(p_rank, p_size, state, prev_full_state);
        if (p_rank == 0 && VERBOSE) {
            for (int i = 0; i < SIZE; i++) {
                printf("%u", prev_full_state[i]);
            }
            printf("\n");
        }
    }

    int iter = 0;
    int converged = 0;
    while (!converged && iter < MAX_ITERS) {
        // share ghost cells
        if (p_size > 1) {
            CHECK_RETURN(sent_recv_ghost_cells(p_rank, p_size, 1 /*to next*/, state), "sent/recv ghost cells to next");
            CHECK_RETURN(sent_recv_ghost_cells(p_rank, p_size, 0 /*to next*/, state), "sent/recv ghost cells to prev");
        }
        else if (CYCLIC_WORLD) {
            for (int i = 0; i < begin(N); i++)
                state[i] = state[end(N) - GHOST_OFFSET + i];
            for (int i = end(N); i < N; i++)
                state[i] = state[GHOST_OFFSET + i - end(N)];
        }

        // calculate next state
        int ghost_lo = state[begin(N) - 1];
        int ghost_hi = state[end(N)];
        converged = make_step(ghost_lo, ghost_hi, begin(N), end(N), state);
        CHECK_RETURN(MPI_Allreduce(
                          &converged,
                          &converged,
                          1,
                          MPI_INT,
                          MPI_LAND,
                          MPI_COMM_WORLD),
                  "MPI_Allreduce");
        iter++;

        // check result with sequential version
        if (DEBUG) {
            gather_states(p_rank, p_size, state, full_state);
            if (p_rank == 0) {
                if (VERBOSE) {
                    for (int i = 0; i < SIZE; i++) {
                        printf("%u", full_state[i]);
                    }
                    printf("\n");
                }
                ghost_lo = CYCLIC_WORLD ? prev_full_state[SIZE-1] : 0;
                ghost_hi = CYCLIC_WORLD ? prev_full_state[0] : 0;
                assert(converged == make_step(ghost_lo, ghost_hi, 0, SIZE, prev_full_state));
                for (int i = 0; i < SIZE; i++) {
                    assert(full_state[i] == prev_full_state[i]);
                }
                memcpy(prev_full_state, full_state, SIZE * sizeof(unsigned int));
            }
        }
    }

    free(state);
    if (p_rank == 0 && DEBUG) {
        free(full_state);
        free(prev_full_state);
    }
    return 0;
}

int main(int argc, char **argv) {
    int p_rank, p_size;
    CHECK_RETURN(MPI_Init(&argc, &argv), "MPI_Init");
    CHECK_RETURN(MPI_Comm_rank(MPI_COMM_WORLD, &p_rank), "MPI_Comm_rank");
    CHECK_RETURN(MPI_Comm_size(MPI_COMM_WORLD, &p_size), "MPI_Comm_size");


    double time_elapsed = MPI_Wtime();
    run_wolfram_automaton(p_rank, p_size);
    time_elapsed = MPI_Wtime() - time_elapsed;
    if (p_rank == 0)
        printf("Time: %lf\n", time_elapsed);

    CHECK_RETURN(MPI_Finalize(), "MPI_Finalize");
    return 0;
}
