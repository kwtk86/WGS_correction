# -*- coding: utf-8 -*-
import os.path

import numpy as np
import inspect
import fiona
import warnings
from typing import Callable

__all__ = ['rec_core', 'rec']

x_pi = 3.14159265358979324 * 3000.0 / 180.0
pi = 3.1415926535897932384626  # π
a = 6378245.0  # 长半轴
ee = 0.00669342162296594323  # 扁率


class CoordTrans:
    LLBAND = [75, 60, 45, 30, 15, 0]
    LL2MC = [
        [-0.0015702102444, 111320.7020616939, 1704480524535203, -10338987376042340, 26112667856603880, -35149669176653700,
         26595700718403920, -10725012454188240, 1800819912950474, 82.5],
        [0.0008277824516172526, 111320.7020463578, 647795574.6671607, -4082003173.641316, 10774905663.51142,
         -15171875531.51559, 12053065338.62167, -5124939663.577472, 913311935.9512032, 67.5],
        [0.00337398766765, 111320.7020202162, 4481351.045890365, -23393751.19931662, 79682215.47186455, -115964993.2797253,
         97236711.15602145, -43661946.33752821, 8477230.501135234, 52.5],
        [0.00220636496208, 111320.7020209128, 51751.86112841131, 3796837.749470245, 992013.7397791013, -1221952.21711287,
         1340652.697009075, -620943.6990984312, 144416.9293806241, 37.5],
        [-0.0003441963504368392, 111320.7020576856, 278.2353980772752, 2485758.690035394, 6070.750963243378,
         54821.18345352118, 9540.606633304236, -2710.55326746645, 1405.483844121726, 22.5],
        [-0.0003218135878613132, 111320.7020701615, 0.00369383431289, 823725.6402795718, 0.46104986909093,
         2351.343141331292, 1.58060784298199, 8.77738589078284, 0.37238884252424, 7.45]]
    # 百度墨卡托转回到百度经纬度纠正矩阵
    MCBAND = [12890594.86, 8362377.87, 5591021, 3481989.83, 1678043.12, 0]
    MC2LL = [[1.410526172116255e-8, 0.00000898305509648872, -1.9939833816331, 200.9824383106796, -187.2403703815547,
              91.6087516669843, -23.38765649603339, 2.57121317296198, -0.03801003308653, 17337981.2],
             [-7.435856389565537e-9, 0.000008983055097726239, -0.78625201886289, 96.32687599759846, -1.85204757529826,
              -59.36935905485877, 47.40033549296737, -16.50741931063887, 2.28786674699375, 10260144.86],
             [-3.030883460898826e-8, 0.00000898305509983578, 0.30071316287616, 59.74293618442277, 7.357984074871,
              -25.38371002664745, 13.45380521110908, -3.29883767235584, 0.32710905363475, 6856817.37],
             [-1.981981304930552e-8, 0.000008983055099779535, 0.03278182852591, 40.31678527705744, 0.65659298677277,
              -4.44255534477492, 0.85341911805263, 0.12923347998204, -0.04625736007561, 4482777.06],
             [3.09191371068437e-9, 0.000008983055096812155, 0.00006995724062, 23.10934304144901, -0.00023663490511,
              -0.6321817810242, -0.00663494467273, 0.03430082397953, -0.00466043876332, 2555164.4],
             [2.890871144776878e-9, 0.000008983055095805407, -3.068298e-8, 7.47137025468032, -0.00000353937994,
              -0.02145144861037, -0.00001234426596, 0.00010322952773, -0.00000323890364, 826088.5]]
    def __init__(self):
        pass

    @staticmethod
    def bd09togcj02(bd_lon, bd_lat):
        """
        百度坐标系(BD09)转火星坐标系(GCJ02)
        :param bd_lat:百度坐标纬度
        :param bd_lon:百度坐标经度
        :return:转换后的坐标列表形式
        """
        x = bd_lon - 0.0065
        y = bd_lat - 0.006
        z = np.sqrt(x * x + y * y) - 0.00002 * np.sin(y * x_pi)
        theta = np.atan2(y, x) - 0.000003 * np.cos(x * x_pi)
        gg_lng = z * np.cos(theta)
        gg_lat = z * np.sin(theta)
        return gg_lng, gg_lat

    @staticmethod
    def gcj02towgs84(lng, lat):
        """
        GCJ02(火星坐标系)转GPS84
        :param lng:火星坐标系的经度
        :param lat:火星坐标系纬度
        :return:
        """
        # if CoordTrans.out_of_china(lng, lat):
        #     return lng, lat
        dlat = CoordTrans.transformlat(lng - 105.0, lat - 35.0)
        dlng = CoordTrans.transformlng(lng - 105.0, lat - 35.0)
        radlat = lat / 180.0 * pi
        magic = np.sin(radlat)
        magic = 1 - ee * magic * magic
        sqrtmagic = np.sqrt(magic)
        dlat = (dlat * 180.0) / ((a * (1 - ee)) / (magic * sqrtmagic) * pi)
        dlng = (dlng * 180.0) / (a / sqrtmagic * np.cos(radlat) * pi)
        mglat = lat + dlat
        mglng = lng + dlng
        return lng * 2 - mglng, lat * 2 - mglat

    @staticmethod
    def bd09towgs84(lng, lat):
        lng, lat = CoordTrans.bd09togcj02(lng, lat)
        return CoordTrans.gcj02towgs84(lng, lat)

    @staticmethod
    def transformlat(lng, lat):
        ret = -100.0 + 2.0 * lng + 3.0 * lat + 0.2 * lat * lat + \
            0.1 * lng * lat + 0.2 * np.sqrt(np.fabs(lng))
        ret += (20.0 * np.sin(6.0 * lng * pi) + 20.0 *
                np.sin(2.0 * lng * pi)) * 2.0 / 3.0
        ret += (20.0 * np.sin(lat * pi) + 40.0 *
                np.sin(lat / 3.0 * pi)) * 2.0 / 3.0
        ret += (160.0 * np.sin(lat / 12.0 * pi) + 320 *
                np.sin(lat * pi / 30.0)) * 2.0 / 3.0
        return ret

    @staticmethod
    def transformlng(lng, lat):
        ret = 300.0 + lng + 2.0 * lat + 0.1 * lng * lng + \
            0.1 * lng * lat + 0.1 * np.sqrt(np.fabs(lng))
        ret += (20.0 * np.sin(6.0 * lng * pi) + 20.0 *
                np.sin(2.0 * lng * pi)) * 2.0 / 3.0
        ret += (20.0 * np.sin(lng * pi) + 40.0 *
                np.sin(lng / 3.0 * pi)) * 2.0 / 3.0
        ret += (150.0 * np.sin(lng / 12.0 * pi) + 300.0 *
                np.sin(lng / 30.0 * pi)) * 2.0 / 3.0
        return ret

    @staticmethod
    def out_of_china(lng, lat):
        """
        判断是否在国内，不在国内不做偏移
        :param lng:
        :param lat:
        :return:
        """
        if lng < 72.004 or lng > 137.8347:
            return True
        if lat < 0.8293 or lat > 55.8271:
            return True
        return False

def trans(features: list, output_shp: str, meta: dict, tran_func: Callable):

    out_dir = os.path.split(output_shp)[0]
    os.makedirs(out_dir, exist_ok=True)

    with fiona.open(output_shp, 'w', **meta, encoding='utf-8') as dst:
        for feature in features:
            if feature['geometry'] is None:
                continue
            parts = feature['geometry']['coordinates']
            if isinstance(parts[0], list): # 多部件形式
                new_parts = []
                for part in parts:
                    part_arr = np.array(part)
                    new_part = tran_func(part_arr[:, 0], part_arr[:, 1])
                    new_part_arr = np.stack(new_part, axis=1)
                    part_arr[:, :2] = new_part_arr
                    new_parts.append(part_arr.tolist())
                feature['geometry']['coordinates'] = new_parts
            elif isinstance(parts[0], tuple): # 单部件
                part_arr = np.array(parts)
                new_part = tran_func(part_arr[:, 0], part_arr[:, 1])
                new_part_arr = np.stack(new_part, axis=1)
                part_arr[:, :2] = new_part_arr
                feature['geometry']['coordinates'] = part_arr.tolist()
            elif isinstance(parts[0], (float, int)): # 点
                new_part = tran_func(parts[0], parts[1])
                part_list = list(parts)
                part_list[:2] = new_part
                feature['geometry']['coordinates'] = part_list
            else:
                raise RuntimeError('Unknown data format')
            dst.write(feature)


def check_func(tran_func: Callable) -> None:
    params_count = len(inspect.signature(tran_func).parameters)
    assert params_count == 2, "自定义函数参数数目必须为2"


def rec_core(input_shp:str, output_shp:str, tran_func: Callable) -> None:
    """
    shapefile坐标纠正核函数，支持自定义转换函数传入
    :param input_shp: 输入shapefile
    :param output_shp: 输出shapefile
    :param tran_func: 坐标转换函数
    :return:
    """
    check_func(tran_func)

    with fiona.open(input_shp, 'r') as src:
        meta = src.meta
        if meta['crs'].is_projected:
            warnings.warn('only WGS 1984(geographical) is supported\n仅支持WGS1984地理坐标系，请检查数据坐标系')
        features = list(src)
        trans(features, output_shp, meta, tran_func)



def rec(input_shp:str, output_shp:str, rec_type: str = 'gd') -> None:
    """
    shapefile坐标纠正
    :param input_shp: 输入shapefile
    :param output_shp: 输出shapefile
    :param rec_type: 转换类型， bd 或 gd
    :return:
    """
    if rec_type == 'bd':
        tran_func = CoordTrans.bd09towgs84
    elif rec_type == 'gd':
        tran_func = CoordTrans.gcj02towgs84
    else:
        raise NotImplementedError('Unknown rectification type')
    rec_core(input_shp, output_shp, tran_func)

