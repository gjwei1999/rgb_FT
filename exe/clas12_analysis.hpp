/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id);
size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id);

size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id) {
  // Called once for each thread
  // Make a new chain to process for this thread
  auto chain = std::make_shared<TChain>("clas12");
  // Add every file to the chain
  for (auto in : inputs) chain->Add(in.c_str());

  // Run the function over each thread
  return run(chain, hists, thread_id);
}

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();
  float beam_energy = NAN;
  if (getenv("CLAS12_E") != NULL) beam_energy = atof(getenv("CLAS12_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  size_t total_twopi = 0;
  size_t total_singlepip = 0;
  
  size_t prot_in_twopi = 0;
  size_t pip_in_twopi = 0;
  size_t pim_in_twopi = 0;
  size_t electron_in_twopi = 0;
  size_t neutron_in_twopi = 0;
  size_t photon_in_twopi = 0;
  size_t others_in_twopi = 0;
  
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
    //  _hists->FillHists_electron_cuts(event, data->ec_tot_energy(0), data->p(0));
    auto event = std::make_shared<Reaction>(data, beam_energy);

 
    // find the position of electron
    int num;

    num = data->gpart();
    int i, k, j;
    k=0;
    j=0;
    for(i=0;i<num;i++){
	if( (data->ft_pid(i) == ELECTRON) && !std::isnan(data->ft_hodo_time(i)) ){        
            data->_pos_of_elec = i;
            k++;
        }
	if(data->ft_pid(i) == ELECTRON){
	j++;
	}
    }
    if(k!=1) continue;   //only 1 electron
    if(j!=1) continue;

    event->SetElec();    

    if (data->charge(data->_pos_of_elec) == -1) _hists->FillHists_elec_cuts(data);
//  if (data->charge(0) == -1) _hists->FillHists_electron_cuts(data);
   
   auto cuts = std::make_shared<Cuts>(data);
  
    if (!cuts->ElectronCuts()) continue;
//    if (!cuts->PipCuts()) continue;

    _hists->Fill_SF(data);
    //_hists->Fill_EC(data->ec_tot_energy(0), data->p(0));

//    // If we pass electron cuts the event is processed
    total++;
      
    // Make a reaction class from the data given
    auto dt = std::make_shared<Delta_T>(data);
       
    // For each particle in the event
    //  if (data->gpart() > 1) continue;
    for (int part = 0; part < data->gpart(); part++) {
      dt->dt_calc(part);
      //
      _hists->Fill_MomVsBeta(data, part);
      _hists->Fill_deltat_pi(data, dt, part);
      _hists->Fill_deltat_prot(data, dt, part);
      //_hists->Fill_deltat_positive(data, dt, event, part);
      _hists->FillHists_elec_with_cuts(data);

      // Check particle ID's and fill the reaction class

      if (cuts->IsProton(part)) {
        event->SetProton(part);
      } else if (cuts->IsPip(part)) {
        event->SetPip(part);
      } else if (cuts->IsPim(part)) {
        event->SetPim(part);
      } else {
        event->SetOther(part);
      }
       event->epsilont();
    }
    // Check the reaction class what kind of even it is and fill the appropriate histograms
    _hists->Fill_histSevenD(event);
    _hists->Fill_WvsQ2(event);
    _hists->Fill_x_mu(event);
    
   // std::cout<<"current event"<<current_event<<std::endl;

    if (event->TwoPion()) {
     // std::cout<<"has twopion reaction"<<std::endl;
      

        for(i=0;i<num;i++){
        if(data->ft_pid(i) ==ELECTRON){
            electron_in_twopi++;
        }
        else if(data->ft_pid(i) ==PROTON){
            prot_in_twopi++;
        }
        else if(data->ft_pid(i) ==PIP){
            pip_in_twopi++;
        }
        else if(data->ft_pid(i) ==PIM){
            pim_in_twopi++;
        }
        else if(data->ft_pid(i) ==NEUTRON){
            neutron_in_twopi++;
        }
        else if(data->ft_pid(i) ==PHOTON){
            photon_in_twopi++;
        }
        else {
            others_in_twopi++;
        }
    }
      
      //if( abs(event->MM2()) >0.05) {continue;}  
      _hists->Fill_WvsQ2_twoPi(event);
      total_twopi++;
    }
    
    if (event->SinglePip()) {
      //std::cout<<"has single pion reaction"<<std::endl;
      _hists->Fill_WvsQ2_singlePi(event);
     total_singlepip++;
    }
    if (event->NeutronPip()) _hists->Fill_WvsQ2_Npip(event);
 


    //
    //
    //modified by Jiawei
    //if (current_event % 100000 == 0) {
    //   std::cout<<current_event<<std::endl;
    //}
    //
    //
    //


}
    std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;

    std::cout << "Twopion Percent = " << 100.0 * total_twopi / num_of_events << std::endl;
    std::cout << "Singlepip Percent = " << 100.0 * total_singlepip / num_of_events << std::endl;

    std::cout << "total number of events = " << num_of_events << std::endl;
    std::cout << "number of twopion events = " << total_twopi << std::endl;
    std::cout << "number of electron in twopion = " << electron_in_twopi << std::endl;
    std::cout << "number of proton in twopion = " << prot_in_twopi << std::endl;
    std::cout << "number of pi+ in twopion = " << pip_in_twopi << std::endl;
    std::cout << "number of pi- in twopion = " << pim_in_twopi << std::endl;
    std::cout << "number of neutron in twopion = " << neutron_in_twopi << std::endl;
    std::cout << "number of photon in twopion = " << photon_in_twopi << std::endl;
    std::cout << "number of others in twopion = " << others_in_twopi << std::endl;
   
    // Return the total number of events
  return num_of_events;
}
#endif
