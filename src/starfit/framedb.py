"""
Describe purpose of this script here

Created: 2/4/25
"""
from contextlib import closing
from os import rename
from pathlib import Path
from sqlite3 import connect
from typing import Callable, Any


def update_db(*,dbname:str,old_fields:dict[str,str],new_fields:dict[str,tuple[str,Callable[[str,dict[str,Any]],Any]]],
              old_tablename:str="frames",new_tablename:str="frames",
              primary_key="framenum"):
    """

    :param old_fields: Dictionary of old fields. Key is name, value is data type
    :param new_fields: Dictionary of new fields. Key is new field name, value is one of the following:
                         * If str: old field name
                         * If callable: Must be a callable equivalent to the following:
                         def f_newvalue(new_fieldname:str,oldrow:dict[str,Any])->Any:
                             '''
                             :param new_fieldname: name of field we are working on
                             :param oldrow: dictionary of old row. Key is old field name, value is field value
                             :return: value for field named new_fieldname for this row in new table
                             '''

    """
    new_table_sql = (f"create table if not exists {new_tablename} ("
                     f"{','.join([k + ' ' + t for k, (t,f) in new_fields.items()])} "
                     f", primary key ({primary_key}))")
    select_sql = (
        f"select {','.join(old_fields.keys())} "
        f"from {old_tablename}")
    insert_sql = (
        f"insert into {new_tablename} ({','.join([new_fieldname for new_fieldname, (t,old_fieldname) in new_fields.items()])}) "
        f"values ({','.join(['?'] * len(new_fields))})")
    old_dbname=dbname.replace(".sqlite","_old.sqlite")
    # rename(dbname,old_dbname)
    with closing(connect(old_dbname)) as old_conn, closing(connect(dbname)) as new_conn:
        with new_conn:
            new_conn.execute(new_table_sql)
        with new_conn:
            with closing(old_conn.cursor()) as old_cur:
                for row in old_cur.execute(select_sql):
                    rowdict={k:v for k,v in zip(old_fields.keys(),row)}
                    this_new_fields={}
                    for new_fieldname,(t,f_newfield) in new_fields.items():
                        new_value=f_newfield(new_fieldname,rowdict)
                        this_new_fields[new_fieldname]=new_value
                    slice_fields=slice(0,len(this_new_fields))
                    this_new_fields={k:v for k,v in list(this_new_fields.items())[slice_fields]}
                    insert_sql=(f"insert or replace into {new_tablename}"
                                f"({','.join([k for k,v in list(this_new_fields.items())])}) "
                                f" values ({','.join(['?' for k,v in list(this_new_fields.items())])})")
                    new_values=list(this_new_fields.values())
                    new_conn.execute(insert_sql,new_values)

def update_db1(casename:str,width:int,height:int,right_num:float):
    dbname=Path(f"data/db/frame_index_{casename}.sqlite")
    old_frames={"framenum":"integer not null",
               "timestamp":"datetime default CURRENT_TIMESTAMP",
               "nstars":"integer",
               "rmsdiff":"real",
               "et":"real",
               "et_sig":"real",
               "et_source":"integer",
               "lat_c":"real",
               "lat_c_sig":"real",
               "lat_c_source":"integer",
               "lon_c":"real",
               "lon_c_sig":"real",
               "lon_c_source":"integer",
               "angle":"real",
               "angle_sig":"real",
               "angle_source":"integer",
               "clock":"real",
               "clock_sig":"real",
               "clock_source":"integer",
               "right_denom":"real",
               "right_denom_sig":"real",
               "right_denom_source":"integer"}
    new_frames={}
    for k,v in old_frames.items():
        new_key=k.replace("_c","")
        new_frames[new_key]=(v,lambda fieldname,row,k=k:row[k])
    # New fields: width and height
    new_frames["width"]=("integer",lambda fieldname,row:width)
    new_frames["height"]=("integer",lambda fieldname,row:height)
    # New field - right numerator, previously an implicit 4.0.
    new_frames["right_num"]=("real",lambda fieldname,row:right_num)
    # Modified field: change right numerator from 4 to match new value
    new_frames["right_denom"]=("real",lambda fieldname,row:row["right_denom"]*(right_num/4) if row["right_denom"] is not None else None)
    new_frames["right_denom_sig"]=("real",lambda fieldname,row:row["right_denom_sig"]*(right_num/4) if row["right_denom_sig"] is not None else None)
    update_db(dbname=str(dbname),old_fields=old_frames,new_fields=new_frames)

def main():
    update_db1("VoyagerUranusHD",width=1280,height=720,right_num=16.0)
    update_db1("VoyagerUranusDetailed",width=720,height=480,right_num=4.0)
    update_db1("VoyagerNeptuneB",width=960,height=720,right_num=4.0)
    update_db1("VoyagerNeptune",width=640,height=480,right_num=4.0)
    update_db1("SuperTrajectory",width=960,height=720,right_num=4.0)
    update_db1("SuperTrajectoryB",width=1280,height=720,right_num=16.0)


if __name__ == "__main__":
    main()
