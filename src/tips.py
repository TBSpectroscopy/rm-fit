#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES), C.P. 160/09
# Universite Libre de Bruxelles
# 50 avenue F. D. Roosevelt, B-1050 Brussels, Belgium
# Phone: +32.2.650.24.18 - E-mail: thibault.bertin@ulb.be - Web: http://www.ulb.ac.be/cpm
#
#-------------------------------------------------------------------------------------------------


def get_tips(tips_file, temperature):

    T = []
    tips = []
    Q = 0.0
    data = False
    with open(tips_file) as f:
        for line in f:
            if data:
                T.append(float(line.split()[0]))
                tips.append(float(line.split()[1]))
                if T[-1] == temperature:
                    Q = tips[-1]
                    return Q
            elif line.startswith("$DATA"):
                data = True

    T_lim = [(0.0, False), (0.0, False)]
    tips_lim = [0.0, 0.0]
    for i in range(0, len(T), 1):
        if T[i] < temperature:
            T_lim[0] = (T[i], True)
            tips_lim[0] = tips[i]
        elif T[i] > temperature:
            T_lim[1] = (T[i], True)
            tips_lim[1] = tips[i]
            break

    if not T_lim[0][1] and not T_lim[1][1]:
        print("ERROR: tips file not covering given temperature")
        return Q

    slope = (tips_lim[1] - tips_lim[0]) / (T_lim[1][0] - T_lim[0][0])
    y0 = tips_lim[0] - (slope * T_lim[0][0])

    Q = y0 + (slope * temperature)
    return Q

                

