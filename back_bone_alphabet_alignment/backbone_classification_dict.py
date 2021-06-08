class BackBone():
    back_bone_dict = {
        "m": "-\u03B3",
        "c": "-\u03C0",
        "l": "-\u03B1",
        "p": "-\u03C1",
        "b": "-\u03B2",
        "t": "\u03B2",
        "y": "\u03C1",
        "a": "\u03B1",
        "i": "\u03C0",
        "g": "\u03B3",
        "n": "-\u03B3"+"i",
        "d": "-\u03C0"+"i",
        "q": "-\u03B1"+"i",
        "r": "-\u03C1"+"i",
        "f": "-\u03B2"+"i",
        "h": "\u03B2"+"i",
        "w": "\u03C1"+"i",
        "k": "\u03B1"+"i",
        "s": "\u03C0"+"i",
        "v": "\u03B3"+"i",
        }
    scoring_matrix = {
        "-γ": {
            "-γ": 7, "-π": 4, "-α": 0, "-ρ": 0, "-β": -1,  "β": 0, "ρ": 0, 
            "α": -1, "π": 3, "γ": 3, "-γi": 4, "-πi": 0, "-αi": 0, "-ρi": 0, 
            "-βi": 0, "βi": 0, "ρi": -2, "αi": 0, "πi": 2, "γi":0
            },
        "-π": {
            "-γ": 4, "-π": 6, "-α": 3, "-ρ": 1, "-β": 0,  "β": 0, "ρ": -1, 
            "α": -1, "π": 1, "γ": 6, "-γi": 8, "-πi": 4, "-αi": 0, "-ρi": -1, 
            "-βi": 3, "βi": 3, "ρi": 0, "αi": 0, "πi": -1, "γi":-1
            },
        "-α": {
            "-γ": 0, "-π": 3, "-α": 4, "-ρ": 2, "-β": 0,  "β": 0, "ρ": -1, 
            "α": -1, "π": 1, "γ": 6, "-γi": 8, "-πi": 4, "-αi": 0, "-ρi": -1, 
            "-βi": 3, "βi": 3, "ρi": 2, "αi": 2, "πi": 2, "γi":1
            }
    }


# From protCAD:
# Backbone Classification Type: -γ -π -α -ρ -β β ρ α π γ -γi -πi -αi -ρi -βi βi ρi αi πi γi
# Backbone Classification Key:   m  c  l  p  b t y a i g  n   d   q   r   f  h  w  k  s  v
