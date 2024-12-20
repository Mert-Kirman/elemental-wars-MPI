# Import the MPI module from the mpi4py library
from mpi4py import MPI

def print_grid(grid):
    '''
    Function that visualizes the current grid
    '''
    if grid is None:
        return

    for i in grid:
        print(' '.join(i))
    print()

if __name__ == "__main__":

    comm = MPI.COMM_WORLD

    # Get the rank (unique ID) of the current process
    rank = MPI.COMM_WORLD.Get_rank()

    # Get the total number of processes in the communicator
    rank_count = MPI.COMM_WORLD.Get_size()

    # Define constants to be used for tags while send, receive
    GRID_TAG = 0
    EARTH_UNIT_COORDINATES_TAG = 1
    FIRE_UNIT_COORDINATES_TAG = 2
    WATER_UNIT_COORDINATES_TAG = 3
    AIR_UNIT_COORDINATES_TAG = 4
    ROW_BOUNDS_TAG = 5
    EXTRA_UPPER_3_ROW_TAG = 6   # From the perspective of the receiving worker process
    EXTRA_LOWER_3_ROW_TAG = 7

    if rank == 0:
        # Read information from input file
        with open("input1.txt", "r") as file:
            lines = file.readlines()

        info = lines[0].strip().split()
        lines = lines[1:]
        
        # Assign simulation parameters
        N = int(info[0])    # Grid size
        W = int(info[1])    # Number of waves
        T = int(info[2])    # Number of units per faction per wave
        R = int(info[3])    # Number of rounds per wave
        
        grid = [["." for i in range(N)] for j in range(N)]

        # Execute simulation
        earth_all = list()  # 2D list that contains all the earth units currently on the grid
        fire_all = list()   # 2D list that contains all the fire units currently on the grid
        water_all = list()
        air_all = list()
        for i in range(W):
            # Keep current locations of each unit in each fraction
            earth_coordinates = lines[i*5+1][2:].strip().split(", ")
            fire_coordinates = lines[i*5+2][2:].strip().split(", ")
            water_coordinates = lines[i*5+3][2:].strip().split(", ")
            air_coordinates = lines[i*5+4][2:].strip().split(", ")
            for j in range(T):
                e_row = int(earth_coordinates[j][0])
                e_col = int(earth_coordinates[j][2])
                if grid[e_row][e_col] == ".":
                    grid[e_row][e_col] = "E"
                    earth_all.append([e_row,e_col])
                
                f_row = int(fire_coordinates[j][0])
                f_col = int(fire_coordinates[j][2])
                if grid[f_row][f_col] == ".":
                    grid[f_row][f_col] = "F"
                    fire_all.append([f_row,f_col])
                
                w_row = int(water_coordinates[j][0])
                w_col = int(water_coordinates[j][2])
                if grid[w_row][w_col] == ".":
                    grid[w_row][w_col] = "W"
                    water_all.append([w_row,w_col])
                
                a_row = int(air_coordinates[j][0])
                a_col = int(air_coordinates[j][2])
                if grid[a_row][a_col] == ".":
                    grid[a_row][a_col] = "A"
                    air_all.append([a_row,a_col])

            # Execute rounds
            for j in range(R):
                print(f'Round {j}:')

                # Prepare grid and other data to send to each worker process
                for k in range(1, rank_count):
                    # Find partition borders for this worker process
                    upper_row_bound = (k-1) * (N // (rank_count-1))
                    lower_row_bound = k * (N // (rank_count-1)) - 1
                    row_bounds = [upper_row_bound, lower_row_bound] # A worker process works between these row bounds

                    # Find which units from each faction will be in this partition
                    earth_partition = [i for i in earth_all if upper_row_bound <= i[0] <= lower_row_bound]
                    fire_partition = [i for i in fire_all if upper_row_bound <= i[0] <= lower_row_bound]
                    water_partition = [i for i in water_all if upper_row_bound <= i[0] <= lower_row_bound]
                    air_partition = [i for i in air_all if upper_row_bound <= i[0] <= lower_row_bound]

                    # Send grid and unit coordinates to worker process
                    comm.send(grid[upper_row_bound:lower_row_bound + 1], dest=k, tag=GRID_TAG)
                    comm.send(earth_partition, dest=k, tag=EARTH_UNIT_COORDINATES_TAG)
                    comm.send(fire_partition, dest=k, tag=FIRE_UNIT_COORDINATES_TAG)
                    comm.send(water_partition, dest=k, tag=WATER_UNIT_COORDINATES_TAG)
                    comm.send(air_partition, dest=k, tag=AIR_UNIT_COORDINATES_TAG)
                    comm.send(row_bounds, dest=k, tag=ROW_BOUNDS_TAG)

                # Receive results from the workers and generate the final grid state
                # TODO
                break


            # print_grid(f'Rank {rank}:\n{print_grid(grid)}')
            break

    else:
        # Receive grid partitions and unit data from the manager
        grid = comm.recv(source=0, tag=GRID_TAG)
        earth_partition = comm.recv(source=0, tag=EARTH_UNIT_COORDINATES_TAG)         
        fire_partition = comm.recv(source=0, tag=FIRE_UNIT_COORDINATES_TAG)         
        water_partition = comm.recv(source=0, tag=WATER_UNIT_COORDINATES_TAG)         
        air_partition = comm.recv(source=0, tag=AIR_UNIT_COORDINATES_TAG)
        row_bounds = comm.recv(source=0, tag=ROW_BOUNDS_TAG)

        upper_row_bound = row_bounds[0]
        lower_row_bound = row_bounds[1]

        # Boundary communication using even-odd communication scheme to avoid deadlocks
        top_boundary = None
        bottom_boundary = None
        if rank > 1:  # Has a top neighbor
            if rank % 2 == 0:  # Even rank: Receive first, then send
                top_boundary = comm.recv(source=rank-1, tag=EXTRA_UPPER_3_ROW_TAG)
                comm.send(grid[:3], dest=rank-1, tag=EXTRA_LOWER_3_ROW_TAG)
            else:  # Odd rank: Send first, then receive
                comm.send(grid[:3], dest=rank-1, tag=EXTRA_LOWER_3_ROW_TAG)
                top_boundary = comm.recv(source=rank-1, tag=EXTRA_UPPER_3_ROW_TAG)

        if rank < rank_count - 1:  # Has a bottom neighbor
            if rank % 2 == 0:  # Even rank: Receive first, then send
                bottom_boundary = comm.recv(source=rank+1, tag=EXTRA_LOWER_3_ROW_TAG)
                comm.send(grid[-3:], dest=rank+1, tag=EXTRA_UPPER_3_ROW_TAG)
            else:  # Odd rank: Send first, then receive
                comm.send(grid[-3:], dest=rank+1, tag=EXTRA_UPPER_3_ROW_TAG)
                bottom_boundary = comm.recv(source=rank+1, tag=EXTRA_LOWER_3_ROW_TAG)

        print(f'Rank={rank}\n')
        print_grid(bottom_boundary)

