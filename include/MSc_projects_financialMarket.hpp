#pragma once

#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <map>

namespace MSC_PROJECTS
{
    enum OrderType { bid , ask };
    /**
     * \brief An order placed in the market
     * \details An order placed on the market, may be to buy/sell. Should include a time stamp, type, price and volume
     */
    struct Order
    {
        /**
         * @details order type
         */
        OrderType oType;
        /**
         * @details time stamp, using integers for ease
         */
        int day;
        /**
         * @details time stamp, using integers for ease
         */
        int period;
        /**
         * @details price, in integer tic levels
         */
        int price;
        /**
         * @details volume of assets requested
         */
        int vol;
        /**
         * @details agent submitting the order
         */
        int agentID;
        /**
         * @brief output an order to a stream
         * @details outputs the data in the format:=
         * ~~~
         *  [BID/ASK] t=(day):(period) p=(price) v=(vol) agentID=(agentID)
         * ~~~
         */
        friend std::ostream& operator<<(std::ostream&,const Order&); 
        
        /**
         * @brief rank orders in the book
         * @details this comparison allows order to be ranked in the book. If the order is an ask(bid), low(high) prices will be ranked first, and further ranked according the day and period submitted so that first come are first served.
         */
        friend bool operator<(const Order& a,const Order&b);
    };
    
    
    struct Agent
    {
        /**
         * @details unique agent ID
         */
        int AgentIndex;
        /**
         * @details number of stocks held by the agent
         */
        int delta;
        /**
         * @details amount of stocks available
         */
        int availableDelta;
        /**
         * @details amount of cash held in the account
         */
        double cash;
        /**
         * @details amount of cash available
         */
        double availableCash;
        /**
         * @details fitness relative to all other agents
         */
        double fitness;
        /**
         * @details agents estimated return
         */
        double r;
        /**
         * @details the identity of the agent that we are following
         */
        int outgoingLink;
        /**
         * @details the identities of the agents that are following us
         */
        std::set<int> incomingLinks;
        /**
         * @details current list of market orders submitted by the agent
         */
        std::set<Order> marketOrders;
        
        /**
         * @details get the current wealth of the agent, given the stock price \f$ S \f$
         * @param St current stock price
         * @return wealth \f$ W_t = B_t + \Delta S_t \f$
         */
        double wealth(double p) const ;
        
        friend std::ostream& operator<<(std::ostream&,const Agent&);
    };
    
    struct Market
    {
        // store the current market price
        double marketPrice;
        
        // omega is used to determine the rate at which new links can be made
        double omega;
        
        // random number generator
        std::mt19937_64 rng;
   
        // store a vector of agents
        std::vector<Agent> marketAgents;
        
        // store the limit order book for ask offers
        std::set<Order> limitOrderBookAsks;
        // store the limit order book for bid offers
        std::set<Order> limitOrderBookBids;
        
        // setup a market with n agents, each holding delta assets and cash
        Market(int n,int delta,double cash):marketAgents(n){
            rng.seed(rng.default_seed);
            for(unsigned int i=0;i<marketAgents.size();i++)
            {
                marketAgents[i].AgentIndex=i;
                marketAgents[i].delta = delta;
                marketAgents[i].availableDelta = delta;
                marketAgents[i].cash = cash;
                marketAgents[i].availableCash = cash;
                marketAgents[i].outgoingLink=i;
            }
        }
        
        // generate an order for a given agent, in a particular time period
        Order strategy(int trader,int period,int tau);
        
        // submit an order to the market
        // - check if trader has sufficient funds/assets
        // - then execute as limit or market order
        int submitOrder(Order o); 
        
        // cancel an order from traders book
        int cancelOrder(Order o); 
        
        // match two orders together
        int matchOrder(const Order &Osubmit,const Order &Omatched);
        
        // execute an order with given volume and price
        int executeOrder(const Order &O,int v,int p); 
        
        // cancel all orders from both trader's and market's books if older than time "s"
        int deleteOrdersMarket(int s); 
        
        // run a simulation with 
        void runSimulation(int totalPeriods, int totalIntraPeriods,double sigma,std::ostream &output=std::cout);
        
        // print full market information
        void printFullMarketData(std::ostream &output=std::cout) const;
        
        // get size of market
        int size() const {return marketAgents.size();}
    };
    
}



