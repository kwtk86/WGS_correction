# -*- coding: utf-8 -*-
import os.path
import numpy as np
import inspect
import fiona
import warnings
from typing import Callable

__all__ = ['correct_core', 'correct']

x_pi = 3.14159265358979324 * 3000.0 / 180.0
pi = 3.1415926535897932384626  # π
a = 6378245.0  # 长半轴
ee = 0.00669342162296594323  # 扁率

class CoordTrans:
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


def correct_core(input_shp:str, output_shp:str, tran_func: Callable) -> None:
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



def correct(input_shp:str, output_shp:str, corr_type: str = 'gd') -> None:
    """
    shapefile坐标纠正
    :param input_shp: 输入shapefile
    :param output_shp: 输出shapefile
    :param corr_type: 转换类型， bd 或 gd
    :return:
    """
    if corr_type == 'bd':
        tran_func = CoordTrans.bd09towgs84
    elif corr_type == 'gd':
        tran_func = CoordTrans.gcj02towgs84
    else:
        raise NotImplementedError('Unknown correction type')
    correct_core(input_shp, output_shp, tran_func)

