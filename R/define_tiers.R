# Specification of tiers and lockdown effect

tier2 = list(res = 0.55, wor = -0.36, gro = -1.41, ret = -1.99, tra = -1.32)
tier3 = list(res = 1.44, wor = -2.72, gro = -1.98, ret = -8.83, tra = -4.62)

lockdown0 = list(res = 0, wor = 0, gro = 0, ret = 0, tra = 0)

# Wales lockdown
#                 variable         V1
# 1:  grocery_and_pharmacy -19.626778
# 2:                 parks -44.768158
# 3:           residential   7.596776
# 4: retail_and_recreation -41.182127
# 5:      transit_stations -20.958923
# 6:            workplaces -21.720779
lockdownW = list(res = 7.59, wor = -21.72, gro = -19.62, ret = -41.18, tra = -20.96)

# NI lockdown
# IMPACT OF LOCKDOWN FOR NORTHERN IRELAND:
#                 variable          V1
# 1:  grocery_and_pharmacy  -0.7792208
# 2:                 parks   7.5283163
# 3:           residential   4.7287157
# 4: retail_and_recreation -13.8701299
# 5:      transit_stations  -9.7984127
# 6:            workplaces -14.4220779
lockdownN = list(res = 4.73, wor = -14.42, gro = -0.78, ret = -13.87, tra = -9.79)

# Scotland restrictions # OLD
# IMPACT OF LOCKDOWN FOR SCOTLAND:
#                 variable         V1
# 1:  grocery_and_pharmacy  0.5931034
# 2:           residential  2.0497126
# 3: retail_and_recreation -7.0221675
# 4:      transit_stations -6.6267492
# 5:            workplaces -8.2720801
lockdownS = list(res = 2.05, wor = -8.27, gro = 0.59, ret = -7.02, tra = -6.63)
