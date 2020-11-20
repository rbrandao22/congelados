#include <algorithm>
#include <cmath>
#include <iostream>
#include <pqxx/pqxx>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

#include "GenArrays.hpp"


// constructor fills products table
GenArrays::GenArrays(const std::vector<std::string> dates_, \
            const std::valarray<double> bins_, const double pop_thres)
{
    pqxx::connection C("dbname = aviacao user = postgres password = passwd"\
            " hostaddr = 127.0.0.1 port = 5432");
    if (C.is_open()) {
        std::cout << "Opened database successfully: " << C.dbname() << \
            std::endl;
    } else {
        std::cout << "Can't open database" << std::endl;
        throw std::runtime_error("aborting");
    }
    pqxx::nontransaction N(C);

    dates = dates_;
    // loop dates
    unsigned counter_dates = 0;
    for (auto& date : dates) {
        // set bin to jan 2020 nominal values
        std::valarray<double> bins = bins_;
        std::string ipca_query = "SELECT index FROM ipca WHERE date = '" + \
                date + "';";
        pqxx::result R_ipca(N.exec(ipca_query));
        auto c = R_ipca.begin();
        bins *= c[0].as<double>() / 100;
    
        // get mkts from given period
        std::string mkts_query = "SELECT DISTINCT origem, destino FROM csv_" + \
                date + ";";
        pqxx::result R_mkts(N.exec(mkts_query));
        for (auto c2 = R_mkts.begin(); c2 != R_mkts.end(); ++c2) {
            std::string pop_query = "SELECT pop_origem, pop_destino FROM " \
                    "mktsinfo WHERE mercado = '" + c2[0].as<std::string>() + \
                    "-" + c2[1].as<std::string>() + "';";
            pqxx::result R_pop(N.exec(pop_query));
            if (R_pop.size() == 0)
                continue;
            auto c3 = R_pop.begin();
            if (c3[0].as<double>() > pop_thres && c3[1].as<double>() > \
                    pop_thres) {
                std::string companies_query = "SELECT DISTINCT empresa FROM "\
                        "csv_" + date + " WHERE origem = '" + c2[0].as\
                        <std::string>() + "' AND destino = '" + c2[1].as\
                        <std::string>() + "';";
                pqxx::result R_companies(N.exec(companies_query));
                if (R_companies.size() > 1) {
                    for (auto c4 = R_companies.begin(); c4 != \
                            R_companies.end(); ++c4) {
                        for (double i = 0; i < bins.size() - 1; ++i) {
                            products.push_back(std::make_tuple\
					       (c2[0].as<std::string>(),\
						c2[1].as<std::string>(),\
						c4[0].as<std::string>(),\
						bins[i], bins[i+1],\
						R_companies.size(),\
						counter_dates));
                        }
                    }
                }
            }
        }

        ++counter_dates;
    }
    std::cout << "GenArrays constructor finished successfully, products " \
            "loaded" << std::endl;
}

void GenArrays::gen_instruments()
{
    // allocate memory
    instruments.resize(products.size());

    // connect to database
    pqxx::connection C("dbname = aviacao user = postgres password = passwd"\
            " hostaddr = 127.0.0.1 port = 5432");
    if (C.is_open()) {
        std::cout << "Opened database successfully: " << C.dbname() << \
                std::endl;
    } else {
        std::cout << "Can't open database" << std::endl;
        throw std::runtime_error("aborting");
    }
    
    pqxx::nontransaction N(C);
    unsigned i = 0;    
    while (i < products.size()) {
        std::string query = "SELECT querosene FROM instruments WHERE date = '"\
	        + dates[std::get<6>(products[i])] + "';";
        pqxx::result R(N.exec(query));
	auto c = R.begin();
	instruments[i] = c[0].as<double>();
        ++i;    
    }    
}

void GenArrays::gen_arrays()
{
    // allocate memory and initialize
    s_obs_wg.resize(products.size());
    mkt_id.resize(products.size());
    std::fill(s_obs_wg.data().begin(), s_obs_wg.data().end(), 0.0);
    pop_ave.resize(products.size());
    X.resize(products.size(), 10);
    // Z matrix w/ 1 instrument, just-identified
    Z.resize(products.size(), 10);

    // connect to database
    pqxx::connection C("dbname = aviacao user = postgres password = passwd"\
            " hostaddr = 127.0.0.1 port = 5432");
    if (C.is_open()) {
        std::cout << "Opened database successfully: " << C.dbname() << \
                std::endl;
    } else {
        std::cout << "Can't open database" << std::endl;
        throw std::runtime_error("aborting");
    }
    pqxx::nontransaction N(C);

    for (auto& date: dates) {
        // loop products
        unsigned i = 0;
        unsigned initial_aux_i = 0;
        unsigned aux_mkt_id = 0;
        double mkt_total = 0;
        std::string curr_mkt = " ";
        std::string mkt = "";
        while (i < products.size()) {
            mkt = std::get<0>(products[i]) + std::get<1>(products[i]);
            if (curr_mkt == " " || mkt == curr_mkt) {
                if (curr_mkt == " ")
                    curr_mkt = mkt;
                std::string query = "SELECT assentos FROM csv_" + date + \
                        " WHERE origem = '" + std::get<0>(products[i]) + \
                        "' AND destino = '" + std::get<1>(products[i]) + \
                        "' AND empresa = '" + std::get<2>(products[i]) + \
                        "' AND tarifa >= " + std::to_string(std::get<3>(products\
                        [i])) + " AND tarifa < " + std::to_string(std::get<4>\
                        (products[i])) + ";";
                pqxx::result R(N.exec(query));
                for (auto c = R.begin(); c != R.end(); ++c) {
                    s_obs_wg[i] += c[0].as<double>();
                }
                mkt_total += s_obs_wg[i];

                // query for pop_ave and X mkt info
                std::string query2 = "SELECT media_g_pop, dist FROM mktsinfo " \
                        "WHERE mercado = '" + std::get<0>(products[i]) + "-" + \
                        std::get<1>(products[i]) + "';";
                pqxx::result R2(N.exec(query2));
                auto c2 = R2.begin();
                // fill pop_ave
                pop_ave[i] = c2[0].as<double>();

                // fill X and Z
                X(i, 0) = 1.;
		Z(i, 0) = 1.;
                X(i, 1) = ((std::get<3>(products[i]) + std::get<4>(products\
                        [i])) / 2.) /1000;
		Z(i, 1) = instruments[i];
                X(i, 2) = c2[1].as<double>() / 1000;
                Z(i, 2) = c2[1].as<double>() / 1000;		
                X(i, 3) = pow(c2[1].as<double>() / 1000, 2);
		Z(i, 3) = pow(c2[1].as<double>() / 1000, 2);
		if (std::get<2>(products[i]) == "AZU") {
		  X(i, 4) = 0.; Z(i, 4) = 0.;
		  X(i, 5) = 0.; Z(i, 5) = 0.;
		  X(i, 6) = 0.; Z(i, 6) = 0.;
		  X(i, 7) = 0.; Z(i, 7) = 0.;
		} else if (std::get<2>(products[i]) == "TAM") {
		  X(i, 4) = 1.; Z(i, 4) = 1.;
		  X(i, 5) = 0.;	Z(i, 5) = 0.;
		  X(i, 6) = 0.;	Z(i, 6) = 0.;
		  X(i, 7) = 0.;	Z(i, 7) = 0.;
		} else if (std::get<2>(products[i]) == "GLO") {
		  X(i, 4) = 0.; Z(i, 4) = 0.;
		  X(i, 5) = 1.;	Z(i, 5) = 1.;
		  X(i, 6) = 0.;	Z(i, 6) = 0.;
		  X(i, 7) = 0.;	Z(i, 7) = 0.;
		} else if (std::get<2>(products[i]) == "ONE") {
		  X(i, 4) = 0.; Z(i, 4) = 0.;
		  X(i, 5) = 0.;	Z(i, 5) = 0.;
		  X(i, 6) = 1.;	Z(i, 6) = 1.;
		  X(i, 7) = 0.;	Z(i, 7) = 0.;
		} else {
		  X(i, 4) = 0.; Z(i, 4) = 0.;
		  X(i, 5) = 0.;	Z(i, 5) = 0.;
		  X(i, 6) = 0.;	Z(i, 6) = 0.;
		  X(i, 7) = 1.;	Z(i, 7) = 1.;
		}
		if (std::get<6>(products[i]) == 0) {
		  X(i, 8) = 0.; Z(i, 8) = 0.;
		  X(i, 9) = 0.; Z(i, 9) = 0.;
		} else if (std::get<6>(products[i]) == 1) {
		  X(i, 8) = 1.; Z(i, 8) = 1.;
		  X(i, 9) = 0.; Z(i, 9) = 0.;
		} else if (std::get<6>(products[i]) == 2) {
		  X(i, 8) = 0.; Z(i, 8) = 0.;
		  X(i, 9) = 1.; Z(i, 9) = 1.;
		}
		
                ++i;
            } else {
                for (unsigned x = initial_aux_i; x <= i; ++x) {
                    s_obs_wg[x] /= mkt_total;
                    mkt_id[x] = aux_mkt_id;
                }
                std::cout << "Finished s, pop_ave and X calcs for mkt " + \
                        curr_mkt << std::endl;
                curr_mkt = mkt;
                mkt_total = 0;
                initial_aux_i = i;
                ++aux_mkt_id;
            }
            if (i == products.size()) {
                for (unsigned x = initial_aux_i; x < i; ++x) {
                    s_obs_wg[x] /= mkt_total;
                    mkt_id[x] = aux_mkt_id;
                }
		std::cout << "Finished s, pop_ave and X calcs for mkt " + \
                        curr_mkt << std::endl;
	    }
        }
    }
}
