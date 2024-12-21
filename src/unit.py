class Unit:
    def __init__(self, unit_type, row, col):
        self.unit_type = unit_type
        self.row = row
        self.col = col
        self.max_health = None
        self.current_health = None
        self.attack_power = None
        self.healing_rate = None
        self.damage_taken = 0

        if unit_type == 'earth':
            self.max_health = 18
            self.current_health = self.max_health
            self.attack_power = 2
            self.healing_rate = 3
        
        elif unit_type == 'fire':
            self.max_health = 12
            self.current_health = self.max_health
            self.attack_power = 4
            self.healing_rate = 1
        
        elif unit_type == 'water':
            self.max_health = 14
            self.current_health = self.max_health
            self.attack_power = 3
            self.healing_rate = 2
        
        elif unit_type == 'air':
            self.max_health = 10
            self.current_health = self.max_health
            self.attack_power = 2
            self.healing_rate = 2


