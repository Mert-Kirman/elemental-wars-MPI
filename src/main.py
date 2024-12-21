# Import the MPI module from the mpi4py library
from mpi4py import MPI
from unit import Unit

def print_grid(grid):
    '''
    Function that visualizes the current grid
    '''
    if grid is None:
        return

    for row in grid:
        for col in row:
            if col == None:
                print('. ')
            else:
                print(f'{col.unit_type[0]} ')
        print()
    print()


def movement_phase(grid, unit_partition, top_boundary, bottom_boundary, upper_row_bound, lower_row_bound, N):
    '''
    Function that moves air units to an adjacent empty cell based on their
    '''
    movement_list = list()
    for unit in unit_partition:
        if unit.type == 'air':
            movement_directions = [(0, 0), (-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
            max_enemy_count = 0
            initial_enemy_count = 0
            max_enemy_row = unit.row
            max_enemy_col = unit.col
            for movement_direction in movement_directions:
                windrush_row, windrush_col = unit.row + movement_direction[0], unit.col + movement_direction[1]
                # Check if out of original grid borders
                if windrush_row < 0 or windrush_row >= N or windrush_col < 0 or windrush_col >= N:
                    continue

                # Check if there are any other units at the windrush location
                new_position = None
                if windrush_row < upper_row_bound:
                    index = upper_row_bound - windrush_row
                    new_position = top_boundary[-index][windrush_col]
                
                elif windrush_row > lower_row_bound:
                    index = windrush_row - lower_row_bound - 1
                    new_position = bottom_boundary[index][windrush_col]

                else:   # Inside partition boundaries
                    index = windrush_row - upper_row_bound
                    new_position = grid[index][windrush_col]

                if new_position != None and new_position.unit_type != 'air': # This air unit can move here
                    continue

                # Find out how many enemies can be attacked at the current windrush position
                attack_directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1), 
                                    (-2, 0), (2, 0), (0, -2), (0, 2), (-2, -2), (-2, 2), (2, -2), (2, 2)]  # Directions an air unit can attack
                local_enemy_count = 0
                for attack_direction in attack_directions:
                    target_row, target_col = windrush_row + attack_direction[0], windrush_col + attack_direction[1]

                    # Check if out of original grid borders
                    if target_row < 0 or target_row >= N or target_col < 0 or target_col >= N:
                        continue

                    # Check if air unit's long range attack is blocked by another unit
                    if attack_direction[0] in [-2, 2] or attack_direction[1] in [-2, 2]:
                        block_target = None
                        is_blocked_target_row = attack_direction[0] // 2
                        is_blocked_target_col = attack_direction[1] // 2
                        if is_blocked_target_row < upper_row_bound:
                            index = upper_row_bound - is_blocked_target_row
                            block_target = top_boundary[-index][is_blocked_target_col]
                        
                        elif is_blocked_target_row > lower_row_bound:
                            index = is_blocked_target_row - lower_row_bound - 1
                            block_target = bottom_boundary[index][is_blocked_target_col]

                        else:   # Inside partition boundaries
                            index = is_blocked_target_row - upper_row_bound
                            block_target = grid[index][is_blocked_target_col]

                        if block_target != None: # Air unit's long range attack is blocked by another unit
                            continue

                    # Check if out of partition borders
                    target = None
                    if target_row < upper_row_bound:
                        index = upper_row_bound - target_row
                        target = top_boundary[-index][target_col]
                    
                    elif target_row > lower_row_bound:
                        index = target_row - lower_row_bound - 1
                        target = bottom_boundary[index][target_col]

                    else:   # Inside partition boundaries
                        index = target_row - upper_row_bound
                        target = grid[index][target_col]

                    if target != None and target.unit_type != 'air': # Enemy found, attack
                        local_enemy_count += 1

                # Calculate how many enemies there are prior to windrush
                if movement_direction == (0, 0):
                    initial_enemy_count = local_enemy_count
                    max_enemy_count = initial_enemy_count
                    max_enemy_row = windrush_row
                    max_enemy_col = windrush_col
                
                # A better attack position using windrush has been found
                elif max_enemy_count < local_enemy_count:
                    max_enemy_count = local_enemy_count
                    max_enemy_row = windrush_row
                    max_enemy_col = windrush_col
                
                # Multiple positions qualify, choose the one with the lowest row-coordinate / column-coordinate
                elif max_enemy_count == local_enemy_count and local_enemy_count != initial_enemy_count:
                    if windrush_row < max_enemy_row:
                        max_enemy_row = windrush_row
                        max_enemy_col = windrush_col
                    elif windrush_row == max_enemy_row and windrush_col < max_enemy_col:
                        max_enemy_row = windrush_row
                        max_enemy_col = windrush_col

            movement_list.append(('move', unit, (max_enemy_row, max_enemy_col)))

    return movement_list


def action_phase(grid, unit_partition, top_boundary, bottom_boundary, upper_row_bound, lower_row_bound, N):
    '''
    Function that decide whether a unit will attack or skip based on their health and strategic considerations
    '''
    action_list = list()
    for unit in unit_partition:
        # Units will skip attacking if their health falls below 50%
        if unit.current_health < unit.max_health // 2:
            action_list.append(('heal', unit))
            continue

        enemy_found = False

        # Attack or skip
        if unit.unit_type == 'earth':
            attack_directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]  # Directions an earth unit can attack
            for attack_direction in attack_directions:
                target_row, target_col = unit.row + attack_direction[0], unit.col + attack_direction[1]
                # Check if out of original grid borders
                if target_row < 0 or target_row >= N or target_col < 0 or target_col >= N:
                    continue
                
                # Check if out of partition borders
                target = None
                if target_row < upper_row_bound:
                    index = upper_row_bound - target_row
                    target = top_boundary[-index][target_col]
                
                elif target_row > lower_row_bound:
                    index = target_row - lower_row_bound - 1
                    target = bottom_boundary[index][target_col]

                else:   # Inside partition boundaries
                    index = target_row - upper_row_bound
                    target = grid[index][target_col]

                if target != None and target.unit_type != 'earth': # Enemy found, attack
                    enemy_found = True
                    action_list.append(('attack', unit, (target_row, target_col)))

                if not enemy_found: # Skip (heal the unit)
                    action_list.append(('heal', unit))

        elif unit.unit_type == 'fire':
            attack_directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]  # Directions a fire unit can attack
            for attack_direction in attack_directions:
                target_row, target_col = unit.row + attack_direction[0], unit.col + attack_direction[1]
                # Check if out of original grid borders
                if target_row < 0 or target_row >= N or target_col < 0 or target_col >= N:
                    continue
                
                # Check if out of partition borders
                target = None
                if target_row < upper_row_bound:
                    index = upper_row_bound - target_row
                    target = top_boundary[-index][target_col]
                
                elif target_row > lower_row_bound:
                    index = target_row - lower_row_bound - 1
                    target = bottom_boundary[index][target_col]

                else:   # Inside partition boundaries
                    index = target_row - upper_row_bound
                    target = grid[index][target_col]

                if target != None and target.unit_type != 'fire': # Enemy found, attack
                    enemy_found = True
                    action_list.append(('attack', unit, (target_row, target_col)))

                if not enemy_found: # Skip (heal the unit)
                    action_list.append(('heal', unit))

        elif unit.unit_type == 'water':
            attack_directions = [(-1, -1), (-1, 1), (1, -1), (1, 1)]  # Directions a water unit can attack
            for attack_direction in attack_directions:
                target_row, target_col = unit.row + attack_direction[0], unit.col + attack_direction[1]
                # Check if out of original grid borders
                if target_row < 0 or target_row >= N or target_col < 0 or target_col >= N:
                    continue
                
                # Check if out of partition borders
                target = None
                if target_row < upper_row_bound:
                    index = upper_row_bound - target_row
                    target = top_boundary[-index][target_col]
                
                elif target_row > lower_row_bound:
                    index = target_row - lower_row_bound - 1
                    target = bottom_boundary[index][target_col]

                else:   # Inside partition boundaries
                    index = target_row - upper_row_bound
                    target = grid[index][target_col]

                if target != None and target.unit_type != 'water': # Enemy found, attack
                    enemy_found = True
                    action_list.append(('attack', unit, (target_row, target_col)))

                if not enemy_found: # Skip (heal the unit)
                    action_list.append(('heal', unit))

        elif unit.unit_type == 'air':
            attack_directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1), 
                                 (-2, 0), (2, 0), (0, -2), (0, 2), (-2, -2), (-2, 2), (2, -2), (2, 2)]  # Directions an air unit can attack

            for attack_direction in attack_directions:
                target_row, target_col = unit.row + attack_direction[0], unit.col + attack_direction[1]
                # Check if out of original grid borders
                if target_row < 0 or target_row >= N or target_col < 0 or target_col >= N:
                    continue
                
                # Check if air unit's long range attack is blocked by another unit
                if attack_direction[0] in [-2, 2] or attack_direction[1] in [-2, 2]:
                    block_target = None
                    is_blocked_target_row = attack_direction[0] // 2
                    is_blocked_target_col = attack_direction[1] // 2
                    if is_blocked_target_row < upper_row_bound:
                        index = upper_row_bound - is_blocked_target_row
                        block_target = top_boundary[-index][is_blocked_target_col]
                    
                    elif is_blocked_target_row > lower_row_bound:
                        index = is_blocked_target_row - lower_row_bound - 1
                        block_target = bottom_boundary[index][is_blocked_target_col]

                    else:   # Inside partition boundaries
                        index = is_blocked_target_row - upper_row_bound
                        block_target = grid[index][is_blocked_target_col]

                    if block_target != None: # Air unit's long range attack is blocked by another unit
                        continue
                    
                # Check if out of partition borders
                target = None
                if target_row < upper_row_bound:
                    index = upper_row_bound - target_row
                    target = top_boundary[-index][target_col]
                
                elif target_row > lower_row_bound:
                    index = target_row - lower_row_bound - 1
                    target = bottom_boundary[index][target_col]

                else:   # Inside partition boundaries
                    index = target_row - upper_row_bound
                    target = grid[index][target_col]

                if target != None and target.unit_type != 'air': # Enemy found, attack
                    enemy_found = True
                    action_list.append(('attack', unit, (target_row, target_col)))

                if not enemy_found: # Skip (heal the unit)
                    action_list.append(('heal', unit))

    return action_list


if __name__ == "__main__":
    comm = MPI.COMM_WORLD

    # Get the rank (unique ID) of the current process
    rank = MPI.COMM_WORLD.Get_rank()

    # Get the total number of processes in the communicator
    rank_count = MPI.COMM_WORLD.Get_size()

    # Define constants to be used for tags while send, receive
    GRID_TAG = 0
    UNIT_TAG = 1
    ROW_BOUNDS_TAG = 2
    EXTRA_UPPER_3_ROW_TAG = 3   # From the perspective of the receiving worker process
    EXTRA_LOWER_3_ROW_TAG = 4
    GRID_SIZE_TAG = 5
    AIR_MOVEMENT_TAG = 6
    ACTION_TAG = 7
    # NEW_EXTRA_UPPER_3_ROW_TAG = 8   # From the perspective of the receiving worker process, new rows after air unit's windrush
    # NEW_EXTRA_LOWER_3_ROW_TAG = 9

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
        
        # grid = [["." for i in range(N)] for j in range(N)]
        grid = [[None for i in range(N)] for j in range(N)]

        # Execute simulation
        units_all = list()  # List that contains all unit objects
        for i in range(W):
            # Keep current locations of each unit in each fraction
            earth_coordinates = lines[i*5+1][2:].strip().split(", ")
            fire_coordinates = lines[i*5+2][2:].strip().split(", ")
            water_coordinates = lines[i*5+3][2:].strip().split(", ")
            air_coordinates = lines[i*5+4][2:].strip().split(", ")
            for j in range(T):
                e_row = int(earth_coordinates[j][0])
                e_col = int(earth_coordinates[j][2])
                if grid[e_row][e_col] == None:
                    new_unit = Unit('earth', e_row, e_col)
                    grid[e_row][e_col] = new_unit
                    units_all.append(new_unit)
                
                f_row = int(fire_coordinates[j][0])
                f_col = int(fire_coordinates[j][2])
                if grid[f_row][f_col] == ".":
                    new_unit = Unit('fire', f_row, f_col)
                    grid[f_row][f_col] = new_unit
                    units_all.append(new_unit)
                
                w_row = int(water_coordinates[j][0])
                w_col = int(water_coordinates[j][2])
                if grid[w_row][w_col] == ".":
                    new_unit = Unit('water', w_row, w_col)
                    grid[w_row][w_col] = new_unit
                    units_all.append(new_unit)
                
                a_row = int(air_coordinates[j][0])
                a_col = int(air_coordinates[j][2])
                if grid[a_row][a_col] == ".":
                    new_unit = Unit('air', a_row, a_col)
                    grid[a_row][a_col] = new_unit
                    units_all.append(new_unit)

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
                    units_partition = [i for i in units_all if upper_row_bound <= i.row <= lower_row_bound]

                    # Send grid and unit coordinates to worker process
                    comm.send(grid[upper_row_bound:lower_row_bound + 1], dest=k, tag=GRID_TAG)
                    comm.send(units_partition, dest=k, tag=UNIT_TAG)
                    comm.send(row_bounds, dest=k, tag=ROW_BOUNDS_TAG)
                    comm.send(N, dest=k, tag=GRID_SIZE_TAG)

                # Receive movement results from the workers and generate the updated grid state after air unit's windrush
                all_movements = list()
                for k in range(1, rank_count):
                    movement_list = comm.recv(source=k, tag=AIR_MOVEMENT_TAG)
                    all_movements.extend(movement_list)

                for movement in all_movements:
                    target_cell = grid[movement[2][0]][movement[2][1]]
                    if target_cell != None: # An air unit exists on this cell in the grid
                        # Combine health
                        tmp_health = target_cell.current_health + movement[1].current_health
                        if tmp_health > target_cell.max_health:
                            tmp_health = target_cell.max_health
                        
                        target_cell.current_health = tmp_health

                        # Combine attack power
                        target_cell.attack_power = target_cell.attack_power + movement[1].attack_power

                        # Delete the remnant of this air unit since it already merged with another unit
                        grid[movement[1].row][movement[1].col] = None

                    else:
                        # Update the air unit object's row and column after windrush
                        grid[movement[1].row][movement[1].col].row = movement[2][0]
                        grid[movement[1].row][movement[1].col].col = movement[2][1]

                        # Move this air unit to its' new location on the grid and delete its' remnant since it moved to a new location
                        grid[movement[2][0]][movement[2][1]] = grid[movement[1].row][movement[1].col]
                        grid[movement[1].row][movement[1].col] = None

                # Update units_all to hold the updated unit objects
                units_all.clear()
                for row in grid:
                    for col in row:
                        if col != None:
                            units_all.append(col)

                # Send new grid and other data info to worker processes after air units complete their windrush
                for k in range(1, rank_count):
                    # Find partition borders for this worker process
                    upper_row_bound = (k-1) * (N // (rank_count-1))
                    lower_row_bound = k * (N // (rank_count-1)) - 1
                    row_bounds = [upper_row_bound, lower_row_bound] # A worker process works between these row bounds

                    # Find which units from each faction will be in this partition
                    units_partition = [i for i in units_all if upper_row_bound <= i.row <= lower_row_bound]

                    # Send grid and unit coordinates to worker process
                    comm.send(grid[upper_row_bound:lower_row_bound + 1], dest=k, tag=GRID_TAG)
                    comm.send(units_partition, dest=k, tag=UNIT_TAG)
                    comm.send(row_bounds, dest=k, tag=ROW_BOUNDS_TAG)
                    comm.send(N, dest=k, tag=GRID_SIZE_TAG)
                
                # Receive action results from the workers and generate the final grid state
                all_actions = list()
                for k in range(1, rank_count):
                    movement_list = comm.recv(source=k, tag=ACTION_TAG)
                    all_actions.extend(movement_list)

                # Calculate damages that will be taken for each unit object
                for action in all_actions:
                    if action[0] == 'attack':
                        attack_power = action[1].attack_power
                        enemy = grid[action[2][0]][action[2][1]]    # Obtain enemy object that will take the damage
                        enemy.damage_taken = attack_power

                # Deal the damages calculated
                for unit in units_all:
                    if unit.unit_type == 'earth':   # Fortification
                        unit.current_health -= unit.damage_taken // 2
                        unit.damage_taken = 0
                    else:
                        unit.current_health -= unit.damage_taken
                        unit.damage_taken = 0

                # Kill units that have current_health <= 0
                units_all.clear()
                units_killed = list()   # Store units that die to apply inferno skills for fire units
                for row in range(N):
                    for col in range(N):
                        if grid[row][col].current_health <= 0: # Kill this unit
                            units_killed.append((row, col))
                            grid[row][col] = None
                        else:
                            units_all.append(grid[row][col])

                # Apply inferno skill for fire units
                # TODO: Increase attack power of inferno units if the unit they attacked is killed

                # TODO: Apply flood skill for water units

                # TODO: Apply heal action

    else:
        # Receive grid partitions and unit data from the manager
        grid = comm.recv(source=0, tag=GRID_TAG)
        unit_partition = comm.recv(source=0, tag=UNIT_TAG)
        row_bounds = comm.recv(source=0, tag=ROW_BOUNDS_TAG)
        N = comm.recv(source=0, tag=GRID_SIZE_TAG)

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

        # Calculate air unit's new positions
        movement_list = movement_phase(grid, unit_partition, top_boundary, bottom_boundary, upper_row_bound, lower_row_bound, N)
        comm.send(movement_list, dest=0, tag=AIR_MOVEMENT_TAG)

         # Receive new grid partitions and unit data from the manager after air's windrush
        grid = comm.recv(source=0, tag=GRID_TAG)
        unit_partition = comm.recv(source=0, tag=UNIT_TAG)

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

        # Calculate actions for each unit in this partition
        action_list = action_phase(grid, unit_partition, top_boundary, bottom_boundary, upper_row_bound, lower_row_bound, N)
        comm.send(action_list, dest=0, tag=ACTION_TAG)
