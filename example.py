import wgs_correction as wcor
wcor.correct(r'route.shp', r'route_1.shp', 'gd')
wcor.correct(r'../data/rec/route_1_1.shp', r'../data/rec/route_3_1.shp', 'gd')
wcor.correct(r'../data/rec/point.shp', r'../data/rec/point_2.shp', 'gd')
wcor.correct(r'../data/rec/county_1.shp', r'../data/rec/county_2.shp', 'gd')