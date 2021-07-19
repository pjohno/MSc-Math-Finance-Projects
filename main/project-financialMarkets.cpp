#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <set>
#include <map>
#include <fstream>


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

/***###################################################################################################################################
***###################################################################################################################################
***###################################################################################################################################
***###################################################################################################################################

  MAIN FUNCTION

***###################################################################################################################################
***###################################################################################################################################
***###################################################################################################################################
***###################################################################################################################################
***###################################################################################################################################
*/

int main()
{
    std::ofstream results("results.csv");
    Market Mfull(100,100,100000);
    Mfull.marketPrice=1000;
    Mfull.runSimulation(200,200,0.05,results);
    Mfull.printFullMarketData();
}

/***###################################################################################################################################
 * 
 * 
 * 
 *  IMPLEMENTATION OF FUNCTIONS
 * 
 * 
 * 
 * 
 */

 
    std::ostream& operator<<(std::ostream& output,const Order& O)
    {
        if(O.oType == bid)
            output << " BID ";
        else if(O.oType == ask)
            output << " ASK ";
        output << " t=" << O.day << ":" << O.period;
        output << " p=" << O.price << " v=" << O.vol << " agentID=" << O.agentID;
        return output;
    }
    
    bool operator<(const Order& a,const Order&b){
        if(a.oType < b.oType)return true;
        else if(a.oType==b.oType)
        {
            if(a.oType == bid && a.price>b.price)
            {
                return true;
            }
            else if(a.oType == ask && a.price<b.price)
            {
                return true;
            }
            else if(a.price==b.price)
            {
                if(a.day<b.day)
                    return true;
                else if(a.day==b.day)
                {
                    if(a.period<b.period)
                    {
                        return true;
                    }
                    else 
                        return false;
                }
                else return false;
                
            }
            else
            {
                return false;
            }
        }
        else return false;
    }
    
    std::ostream& operator<<(std::ostream& output,const Agent& A)
    {
        output << " Agent:= " << A.AgentIndex << std::endl;
        output << " Cash:= " << A.cash << " Assets:= "<< A.delta << std::endl;
        output << " Available Cash:= " << A.availableCash << " Available Assets:= "<< A.availableDelta << std::endl;
        //         output << " Order Book:= \n";
        int runningDelta= A.delta;
        int runningCash=A.cash;
        for(const auto& o : A.marketOrders)
        {
            output << " => " << o << "\n";
            if(o.oType == bid)
                runningCash = runningCash - o.vol*o.price;
            else
                runningDelta = runningDelta - o.vol;
        }
        if(fabs(runningCash - A.availableCash)>1.e-8)
            output << "## ERROR :: CASH NOT ADDING UP " << std::endl;
        if(runningDelta!=A.availableDelta)
            output << "## ERROR :: ASSETS NOT ADDING UP " << std::endl;        
        return output;
    }
    
    void Market::runSimulation(int totalPeriods, int totalIntraPeriods,double sigma,std::ostream &output)
    {
        output << "#period  , #marketPrice" << "\n";          
        // run a bunch of 
        for(int period=1;period<=totalPeriods;period++)
        {
            // calculate the current price of the asset from order books
            if(limitOrderBookAsks.size()>0 && limitOrderBookBids.size()>0 )
            {
                double bestBidPrice = limitOrderBookBids.begin()->price;
                double bestAskPrice = limitOrderBookAsks.begin()->price;
                marketPrice = ( bestAskPrice + bestBidPrice ) / 2.;
            }
            
            // update agent expectations
            for(auto& A : marketAgents)
            {
                std::normal_distribution<double> eps(0.,sigma);
                A.r = eps(rng);
            }
            
            for(int tau = 0;tau<=totalIntraPeriods;tau++)
            {
                // randomly select an agent
                std::uniform_int_distribution<> U(0,size()-1);
                int randomlySelectedAgent = U(rng);
                Order o = strategy(randomlySelectedAgent,period,tau);
                submitOrder(o); 
            }
            
            output << period << " , " << marketPrice << "\n";          
            deleteOrdersMarket(period-2);
            
        }
        output << std::endl;          
            
    }
    
    Order Market::strategy(int agentID,int day, int period)
    {
        Order O;
        Agent& trader = marketAgents[agentID];
        // calculate price
        double pt=marketPrice*exp(trader.r);
        // set time stamp 
        O.day = day;O.period=period;
        // and unique ID
        O.agentID = trader.AgentIndex;
        std::uniform_real_distribution<> U(0.,1.);
        if(pt>marketPrice)
        {
            O.oType = bid;
            O.price = 0;
            double xi1=U(rng),xi2=U(rng);
            int pT = marketPrice*(1-xi1)*(1-xi2) + pt*xi2 + 0.5;
            O.price=std::max(pT,O.price);
            // can only commit to buying given what you've already bud for, so check available cash
            std::uniform_real_distribution<> V(0.,std::max(0.,trader.availableCash)/O.price);
            O.vol = int( V(rng) );
        }
        else
        {
            O.oType = ask;  
            O.price = 2*pt;
            double xi1=U(rng),xi2=U(rng);
            int pT = pt*(1-xi2) + marketPrice*(1+xi1)*xi2 + 0.5;
            O.price=std::min(pT,O.price);
            // can only sell assets you have, so check trader's available assets
            std::uniform_int_distribution<> V(0,std::max(0,trader.availableDelta));
            O.vol = int( V(rng) );
        }
        return O;
        
    }
    
    int Market::submitOrder(Order o)
    {
        // get the trader ID from the order
        Agent& trader = marketAgents[o.agentID];
        // check if order already exists, if so exit now
        auto oo = trader.marketOrders.find(o);
        if(oo!=trader.marketOrders.end())return -1;
       
        // first add order to the agents book, 
        if(o.vol>0)
        {
            // adjust available assets at this stage
            if(o.oType== bid)
            {
                double cash=trader.availableCash - o.vol*o.price;
                // if cash not available to execute transaction then exit
                if(cash>=0.)
                    trader.availableCash = cash;
                else
                    return -1;
            }
            else
            {
                double vol= trader.availableDelta - o.vol;
                // if assets not available then exit
                if(vol>=0)
                    trader.availableDelta = vol;
                else 
                    return -1;
            }
            // OK this is valid order
            // put it in the market orders
            trader.marketOrders.insert(o);
        }
        else
        {
            // return from function if agent selects volume = 0  
            return -1;
        }
        
        // next track the order on the market books
        
        // create a copy of the order to track volumes
        Order Onew(o);
        
        // now look at bid offer first
        if(Onew.oType == bid)
        {
            // iterate through the limit order book of ask offers,
            // to see if any match the bid
            // askOIt will be the current best ask 
            auto askOIt = limitOrderBookAsks.begin();
            
            // go through the book
            while( askOIt != limitOrderBookAsks.end() )
            {
                // get a copy of the current best ask order
                Order askCurrentBest(*askOIt);
                // if there are any ask offers below the bid then try and match them
                if( askCurrentBest.price <= Onew.price )
                {
                    // EXECUTE ORDER
                    int volumeRemaining = matchOrder(Onew,askCurrentBest);
                    //                     std::cout << " VOLUME REMAINING " << volumeRemaining << std::endl;
                    if(volumeRemaining<0)
                    {
                        std::cout << "## error matching agents " << std::endl;
                    }
                    
                    // trade executed
                    int volumeExecuted = Onew.vol - volumeRemaining;
                    if(volumeExecuted>0)marketPrice = askCurrentBest.price;
                    // remove executed trades from ask on book
                    askCurrentBest.vol = ( askCurrentBest.vol - volumeExecuted );
                    // remove the ask order from the book
                    askOIt = limitOrderBookAsks.erase(askOIt);
                    
                    // if volume is unexecuted on the ask side,
                    // add it back in
                    if(askCurrentBest.vol>0)
                        limitOrderBookAsks.insert(askCurrentBest);
                    
                    // now check the bid
                    Onew.vol = volumeRemaining;
                    // EXECUTE ORDER FOR BOTH AGENTS
                    if(Onew.vol==0)
                    {
                        // 
                        cancelOrder(o);
                        return 0;
                    }
                }
                else
                {
                    // this mean best ask is above bid, so break from loop
                    break;
                }
            }
            // if we reach this point, it is because there is still some volume left on the bid unmatched
            // if so, add it to the book
            limitOrderBookBids.insert(Onew);
        }
        else
        {
            // iterate through the limit order book of bid offers,
            // to see if any match the ask
            // bidOIt will be the current best ask 
            auto bidOIt = limitOrderBookBids.begin();
            
            // run through bids
            while( bidOIt != limitOrderBookBids.end() )
            {
                Order bidCurrentBest(*bidOIt);
                // if there are any ask offers below the bid then try and match them
                if( bidCurrentBest.price >= Onew.price )
                {
                    // EXECUTE ORDER
                    int volumeRemaining = matchOrder(Onew,bidCurrentBest);
                    if(volumeRemaining<0)
                    {
                        std::cout << "## error matching agents " << std::endl;
                    }
                    
                    // trade executed
                    int volumeExecuted = Onew.vol - volumeRemaining;
                    if(volumeExecuted>0)marketPrice = bidCurrentBest.price;
                    bidCurrentBest.vol = ( bidOIt->vol - volumeExecuted );
                    
                    //                     std::cout << " ERASE OLD BID " << std::endl;
                    bidOIt = limitOrderBookBids.erase(bidOIt);
                    
                    // if there is any volume left on the current best bid, add it back into market
                    if(bidCurrentBest.vol>0)
                    {
                        limitOrderBookBids.insert(bidCurrentBest);
                    }
                    
                    Onew.vol = volumeRemaining;
                    // RETURN 0 if order executed by market book
                    //                     std::cout << " RETURN 0 if order executed by market book " << std::endl;
                    if(Onew.vol==0){
                        cancelOrder(o);
                        return 0;
                    }
                    //                     std::cout << " MOVE TO NEXT BID IN THE BOOK " << std::endl;
                    
                }
                else
                {
                    // 
                    break;
                }
            }
            //             std::cout << " New ask :: " << Onew << std::endl;
            limitOrderBookAsks.insert(Onew);
        }
        return Onew.vol;
    }
    
    int Market::cancelOrder(Order O)
    {
        auto& A = marketAgents[O.agentID];
        // check if the order is logged in the agents book
        auto oo = A.marketOrders.find(O);
        // if we can find it remove it
        if(oo != A.marketOrders.end())
        {
            // first update available cash/stock to reflect removal of order
            if(O.oType == bid)
            {
                // order is a bid so we return available cash
                A.availableCash = A.availableCash + O.vol * O.price;
            }
            else
            {
                // order is a ask so we return available stock
                A.availableDelta = A.availableDelta + O.vol;
            }
            // then delete it
            A.marketOrders.erase(O);
            // return 0 if order cancelled by agent
            return 0;
        }
            // return -1 if agent doesn't have this order  
        return -1;
    
    }

    int Market::deleteOrdersMarket(int s)
    {
        auto bidOIt = limitOrderBookBids.begin();
        while(bidOIt != limitOrderBookBids.end())
        {
            if(bidOIt->day < s)
            {
                cancelOrder(*bidOIt);
                bidOIt = limitOrderBookBids.erase(bidOIt);
            }
            else
                bidOIt++;
        }
        auto askOIt = limitOrderBookAsks.begin();
        while(askOIt != limitOrderBookAsks.end())
        {
            if(askOIt->day < s)
            {
                cancelOrder(*askOIt);
                askOIt = limitOrderBookAsks.erase(askOIt);
            }
            else
                askOIt++;
        }
        return 0;
    }

    
    int Market::matchOrder(const Order &oExe,const Order &oCou)
    {
        // check bid/ask is opposite
        if(oExe.oType == oCou.oType )
        {
            std::cout << "## can't match 2 bids or 2 asks " << std::endl;
            return -1;
        }
        // try to execute v shares at order price
        int unmatchedVol = executeOrder(oCou,oExe.vol,oCou.price);
        if(unmatchedVol<-1)
        {
            std::cout << "## Agent #" << oCou.agentID;
            std::cout << "## has no matching order." << std::endl;
        }
        // volume executed
        int vExecuted = oExe.vol - unmatchedVol;
        // execute executive order
        executeOrder(oExe,vExecuted,oCou.price);
        return unmatchedVol;
    }
    
    
    int Market::executeOrder(const Order &O,int v,int p)
    {
       // get trader
        auto& trader = marketAgents[O.agentID];
        // execute the given order O at the price p, with a volume of at most v
        // check that this agent has the order in his book
        std::set<Order>::iterator oo = trader.marketOrders.find(O);
        // if valid execute order and delete from the agents book
        if(oo != trader.marketOrders.end())
        {
            // try executing minimum volume from Order and requested v
            int volume = std::min(O.vol,v);
//              std::cout << "## EXECUTE " << O << " with vol " << volume << " at price " << p << std::endl;
            // order type
            if(O.oType == bid)
            {
                // check their is sufficient cash to execute. If not execute maximum amount possible.
                volume = std::min(volume,int(trader.cash/p));
                // order is a bid so calculate transaction cost and take it from cash
                double transactionCost = volume * p;
                trader.cash = trader.cash - transactionCost;
                // add stocks 
                trader.delta = trader.delta + volume;
                
                // will be cancelling this order now, so add back in total order cost, and remove actual transaction
                trader.availableCash = trader.availableCash + O.vol * O.price - transactionCost;
                // add stocks to available
                trader.availableDelta = trader.availableDelta + volume;
            }
            else 
            {
                // check delta is sufficiently large to execute. If not execute maximum amount possible.
                volume = std::min(volume,trader.delta);
                // order is an askso calculate transaction cost and take it from cash
                double transactionCost = volume * p;
                trader.cash = trader.cash + transactionCost;
                trader.delta = trader.delta - volume;
                
                // will be cancelling this order now, so add back in total stocks in order, before removing actual amount sold
                trader.availableCash = trader.availableCash + transactionCost;
                trader.availableDelta = trader.availableDelta + O.vol - volume;
            }
            // Create a new order and adjust the volume to the remaining left over after execution
            Order Onew(O);
            Onew.vol = Onew.vol - volume;
            // now remove it from my order list
            trader.marketOrders.erase(oo);
            
            // if any outstanding volume left, add it back in
            // put it in the market orders
            if(Onew.vol>0){
                trader.marketOrders.insert(Onew);
                // adjust available assets at this stage
                if(Onew.oType== bid)
                    trader.availableCash = trader.availableCash - Onew.vol*Onew.price;
                else
                    trader.availableDelta = trader.availableDelta - Onew.vol;
            }
            // return number of orders left unmatched
            return v-volume;
        }
        // return -1 if agent doesn't have this order  
        return -1;
    }
    
    void Market::printFullMarketData(std::ostream& output) const
    {
        output << "###############\n# MARKET DATA\n###############\n";
        output << " Market Price " << marketPrice << "\n";
        output << "#######\n# MARKET BIDS\n######\n";
        for(auto b : limitOrderBookBids)output << b << "\n";
        output << "#######\n# MARKET ASKS\n######\n";
        for(auto a : limitOrderBookAsks)output << a << "\n";
        output << "#######\n# MARKET AGENTS\n######\n";
        for(auto A : marketAgents)output << A << "\n";

        
    }
