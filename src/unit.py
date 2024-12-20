class Unit:
    def __init__(self, unit_type, row, col):
        self.unit_type = unit_type
        self.row = row
        self.col = col
        self.health = None
        self.attack_power = None
        self.healing_rate = None

        if unit_type == 'earth':
            self.health = 18
            self.attack_power = 2
            self.healing_rate = 3
        
        elif unit_type == 'fire':
            self.health = 12
            self.attack_power = 4
            self.healing_rate = 1
        
        elif unit_type == 'water':
            self.health = 14
            self.attack_power = 3
            self.healing_rate = 2
        
        elif unit_type == 'air':
            self.health = 10
            self.attack_power = 2
            self.healing_rate = 2


