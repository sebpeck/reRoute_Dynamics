import importlib
from reRoute_Dynamics import Physics_Engine as pe
from reRoute_Dynamics import Geography_Tools as gt
from reRoute_Dynamics import Object_Params as op
from reRoute_Dynamics import Trip_Simulator as ts
import multiprocessing
import random
import pandas as pd
import numpy as np
import geopandas as gpd
import gc

importlib.reload(pe)
importlib.reload(ts)

def process_result_list(res_list, ess):
    '''
    Method for processing a list of outputs from ts.run_trip()
    '''
    def subprocess(res, ess):
        res = pd.DataFrame(res)
        bp = res['BP']
        volt = res['v']
        
        milage = res['dx'].sum()/1609.344 # to miles
        time_elapsed = res['dt'].sum() # in seconds
        net_en = (bp/1000*res['dt']).sum()/3600
        SOH_Loss = res['SOH_loss'].sum()
        
        return {"miles":milage,
                "kwh":net_en,
                "tot_time":time_elapsed,
                "v":volt,
                "dt":res['dt'],
                'g':res['grade'],
                'vel':res['v_f'],
                'dx':res['gdx'],
                'geo':res['geometry'],
                'P':res['P'],
                'BP':res['BP'],
                'SOH_loss':res['SOH_loss']}
    
    dat = pd.Series(res_list)
    out = dat.apply(lambda x: subprocess(x, ess))

    return out


def run_shortname_MC(route_shortname, available_routes, seeds, trip_path, traffic_min_max, signal_min_max, wind_speed_min_max, wind_directions_min_max, t_sigs, t_stops, ESS_path, p_aux_min_max, b_factors, a_factors, max_dist, batch_size):
    print("-----------------------{}-------------------".format(route_shortname))

    selected_route_instances = available_routes.get_group(route_shortname).reset_index(drop=True)
    selected_route_instances['results'] = [[]]*len(selected_route_instances)

    
    route_frame = available_routes.get_group(route_shortname)

    route_frame = route_frame.groupby('direction').get_group(1)
    route_frame['bus'] = len(route_frame)*[['./KC_Example_Data/Saved_Objects/Busses/XDE60.txt','./KC_Example_Data/Saved_Objects/Busses/XDE40.txt','./KC_Example_Data/Saved_Objects/Busses/XDE35.txt']]
    
    
    selected_route_instances = route_frame.explode('bus').reset_index()
    del route_frame
    gc.collect()
    selected_route_instances['results'] = [[]]*len(selected_route_instances)


    # Simulation setup
    sim_name = ""
    for k in range(len(selected_route_instances)):
        for index, row in selected_route_instances[k:k+1].iterrows():
            sim_name = '{}_{}_{}_{}'.format(route_shortname, row['shape'], row.direction, row.bus.split('/')[-1].split('.')[0])
            route = gt.load_from_json(row.path)
            busses= []
            trips = []
            esses= []
            for seed in seeds:
                random.seed(int(seed))
                bus = op.load_bus_params(row.bus)
                trip = op.load_trip_params(trip_path)
                trip.m_riders = random.uniform(*row.ridership_min_max)
                trip.traffic = random.uniform(*traffic_min_max)
                trip.chance_sig = random.uniform(*signal_min_max)
                trip.v_wind = random.uniform(*wind_speed_min_max)
                trip.wind_bearing = gt.compass_heading(random.randrange(*wind_directions_min_max))
                trip.seed=int(seed)
                trip.t_sig = random.uniform(*t_sigs)
                trip.t_stop = random.uniform(*t_stops)
                trips.append(trip)
                ESS_path = ESS_path[:-13] + row.bus.split("/")[-1].split('.')[0] + "_ESS.txt"  #Bad code, made to work and get proper associated ESS
                ess = op.load_ESS_params(ESS_path)
                ess.P_aux = random.randrange(*(p_aux_min_max))
                esses.append(ess)
                bus.f_br = random.uniform(*b_factors)
                bus.f_a = random.uniform(*a_factors)
                bus.dmax = random.uniform(*max_dist)
                busses.append(bus)

            all_runs = list(zip([route]*len(trips), trips, busses, esses))
            batched_runs = [all_runs[i:i + batch_size] for i in range(0, len(all_runs), batch_size)]
            print("[{}] : Simulation parameters set.".format(sim_name))

            res_list = []
            for i in range(len(batched_runs)):
                #print(index/len(selected_route_instances[k:k+1])*100, i/len(batched_runs)*100, end= '\r')
                print("[{}] : Batch {} of {}".format(sim_name, i, len(batched_runs)), end='\r')
                batch = batched_runs[i]
                with multiprocessing.Pool(batch_size) as pool:
                    res_list.extend(pool.starmap(ts.simulate_trip, batch))


            proc_res = process_result_list(res_list, ess)
            
            print("[{}] : Results Processed.".format(sim_name))

            out_df = pd.DataFrame(list(proc_res))
            out_df['v'] = out_df['v'].apply(lambda x: x.to_dict())
            out_df['dt'] = out_df['dt'].apply(lambda x: x.to_dict())
            out_df['g'] = out_df['g'].apply(lambda x: x.to_dict())
            out_df['vel'] = out_df['vel'].apply(lambda x: x.to_dict())
            out_df['dx'] = out_df['dx'].apply(lambda x: x.to_dict())
            out_df['geo'] = out_df['geo'].apply(lambda x: x.to_dict())
            out_df['P'] = out_df['P'].apply(lambda x: x.to_dict())
            out_df['SOH_loss'] = out_df['SOH_loss'].apply(lambda x: x.to_dict())
            out_df['95'] = pd.Series(row.results).apply(lambda x: x['kwh']).quantile(.95)
            out_df['97'] = pd.Series(row.results).apply(lambda x: x['kwh']).quantile(.97)
            out_df['99'] = pd.Series(row.results).apply(lambda x: x['kwh']).quantile(.99)
            out_df['99.9'] = pd.Series(row.results).apply(lambda x: x['kwh']).quantile(.999)
            out_df.to_pickle('../../rRD-MISC/MC_Saves_6/{}_{}_{}_{}.pk'.format(route_shortname, row['shape'], row.direction, row.bus.split('/')[-1].split('.')[0]))
            print("[{}] : File Saved.".format(sim_name))
            del out_df
            gc.collect()


    gc.collect()
