from matplotlib import pyplot as plt
from scipy.optimize import minimize
from scipy.special import expit
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Qt5Agg")
import os
import json

import psytrack_learning as psy
from psytrack_learning.getMAP import getMAP
from psytrack_learning.helper.helperFunctions import update_hyper, hyper_to_list
from psytrack_learning.helper.jacHessCheck import compHess, compHess_nolog
from psytrack_learning.helper.invBlkTriDiag import getCredibleInterval
from psytrack_learning.hyperparameter_optimization import evd_lossfun
from psytrack_learning.learning_rules import RewardMax, PredictMax, REINFORCE, REINFORCE_base
from psytrack_learning.simulate_learning import reward_max, predict_max, reinforce, reinforce_base
from psytrack_learning.simulate_learning import simulate_learning

plt.ion()

# Set matplotlib defaults from making files editable and consistent in Illustrator
colors = psy.COLORS
zorder = psy.ZORDER
plt.rcParams['figure.dpi'] = 140
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.facecolor'] = (1,1,1,0)
plt.rcParams['savefig.bbox'] = "tight"
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['axes.labelsize'] = 12

# fit to all AB first trial then save the AIC and BIC

root_dir = r'Z:\HongliWang\Juvi_ASD Deterministic'
strains = ['ChD8']

protocols = ['AB', 'AB-CD']
nSession = 3

dataIndexPath = os.path.join(root_dir, strains[0], 'dataIndex.csv')
dataIndex = pd.read_csv(dataIndexPath)
subjects = dataIndex['Animal'].unique()

# loop through all 6 sessions to get a result
for pp in protocols:
    for nSes in range(1,nSession+1):
        ABList = []
        subject_inuse = []
        for sub in subjects:
            mask = (
                    (dataIndex["Animal"] == sub) &
                    (dataIndex["Protocol"] == pp) &
                    (dataIndex["ProtocolDay"] == nSes)
                    )
            if sum(mask)>0:
                subject_inuse.append(sub)
                filename = dataIndex.BehCSV[mask]
                ABList.append(os.path.join(root_dir, strains[0], 'Analysis',str(sub), filename.iloc[0]))
        # filename)

#%% fit the model, save the AIC and BIC results

#%% REINFORCE learning rule
        AICList = np.zeros((len(ABList),1))
        BICList = np.zeros((len(ABList),1))
        negLLList = np.zeros((len(ABList),1))
        if pp == 'AB':
            sessionLabel = 'AB-'+pp+str(nSes)
        elif pp=='AB-CD':
            sessionLabel = pp+'-CD'+str(nSes)

        for idx,ff in enumerate(ABList):
            #filename1 = r'Z:\HongliWang\Juvi_ASD Deterministic\TSC2\Analysis\361\361_20231119_behaviorDF.csv'
            #filename2 = r'Z:\HongliWang\Juvi_ASD Deterministic\TSC2\Analysis\317\317_20230710_behaviorDF.csv'
            #filename3 = r'Z:\HongliWang\Juvi_ASD Deterministic\TSC2\Analysis\317\317_20230711_behaviorDF.csv'
            #data1 = pd.read_csv(filename1)
            #data2 = pd.read_csv(filename2)
            #data3 = pd.read_csv(filename3)
            #data = pd.concat([data1, data2,data3])
            data = pd.read_csv(ff)
            #data = data1
            subject= subject_inuse[idx]

            # remove nans
            data.replace({"actions": ["NAN","NaN", "nan", "None", ""]}, np.nan, inplace=True)
            data.replace({"schedule": ["NAN", "NaN", "nan", "None", ""]}, np.nan, inplace=True)
            data = data.dropna(subset=["actions"])
            data = data.dropna(subset=['schedule'])
            data.schedule = data.schedule.astype(int)
            # if CD session
            if 'CD' in sessionLabel:
                # remove AB trials
                data = data[np.logical_and(data['schedule']!=1,data['schedule']!=2)].copy()
                data['schedule'] = data['schedule']-2

            s_a = data['schedule'].copy(); s_b = data['schedule'].copy()
            s_a[s_a==2] = 0
            s_b[s_b==1] = 0
            s_b[s_b==2] = -1
            correct = data['reward'].copy()
            correct[np.isnan(correct)] = 0
            correct[correct>0] = 1
            answer = (data['schedule'].copy()-1.5)*2
            y = data['actions'].copy()
            dayLength = np.array([1])

            # determine inputs
            t = np.array(data["trial"])
            prior = ((t[1:] - t[:-1]) == 1).astype(int)
            prior = np.hstack(([0], prior))

            # Calculate previous average tone value
            s_avg = (s_a + s_b)/2
           # s_avg = (s_avg - np.mean(s_avg))/np.std(s_avg)
           # s_avg = np.hstack(([0], s_avg))
           # s_avg = s_avg * prior  # for trials without a valid previous trial, set to 0

            # Calculate previous correct answer
            h = correct[:-1].astype(int)   # map from (0,1) to (-1,1)
            h = np.hstack(([0], h))
            h = h * prior  # for trials without a valid previous trial, set to 0

            # Calculate previous choice
            c = y[:-1].astype(int)   # map from (0,1) to (-1,1)
            c = np.hstack(([0], c))
            c = c * prior  # for trials without a valid previous trial, set to 0
            stick_all = np.insert(np.array((y[0:-1]-0.5)*2),0,0) # do not depend on reward
            stick = np.insert(np.array(correct[0:-1])*np.array((y[0:-1]-0.5)*2),0,0)
            inputs = dict(s_a = np.array(s_a)[:, None],
                          s_b = np.array(s_b)[:, None],
                          cBoth = np.array(s_avg)[:, None],
                          stick = stick[:,None],
                          stick_all=stick_all[:,None],
                          h = np.array(h)[:, None],
                          c = np.array(c)[:, None])


            dat = dict(
                subject = subject_inuse,
                inputs = inputs,
                s_a = np.array(s_a),
                s_b = np.array(s_b),
                correct = np.array(correct),
                answer = np.array(answer),
                y =np.array(y),
                dayLength=dayLength,
            )

            seed = 42
            #mouse_name = "CSHL_003"
            np.random.seed(seed)

            #orig_dat = getMouse(mouse_name, 5)
            orig_dat = dat
            #trim_dat = psy.trim(orig_dat, END=N)
            trim_dat = dat.copy()

            # weights = {"bias": 1, "cBoth": 1}
            # K = np.sum([i for i in weights.values()])
            # X = psy.read_input(trim_dat, weights)
            # answer = trim_dat["answer"]
            # true_sigma  = 2**np.array([-4.0, -5.0])
            # true_alpha  = 2**np.array([-6.0, -8.0])
            # sim_learning_rule = reinforce
            rec_learning_rule = REINFORCE
            #
            # W, y, r, sim_noise = simulate_learning(X=X, answer=answer, sigma=true_sigma, sigma0=1,
            #                                        alpha=true_alpha, learning_rule=sim_learning_rule)
            #
            # dat = {'inputs': trim_dat['inputs'].copy(), 'y': y, 'answer': answer, "correct": r}
            # gen_dat = {"dat": dat, 'true_sigma': true_sigma, 'true_alpha': true_alpha, "weights": weights, "K": K,
            #            'sim_learning_rule': sim_learning_rule, "rec_learning_rule": rec_learning_rule,
            #            "W": W, "sim_noise": sim_noise, "seed": seed}
            # fig = psy.plot_weights(W, weights)
            # plt.show()
            #
            # hyper_guess = {
            #     'alpha': [2**-6] * K,
            #     'sigma': [2**-4] * K,
            #     'sigInit': [2**4] * K,
            #     'sigDay': None,
            # }
            # rec_learning_rule = REINFORCE
            # # Optimizing for both sigma and alpha simultaneously
            # optList = ['sigma', 'alpha']
            #
            # # List of extra arguments used by evd_lossfun in optimization of evidence
            # args = {"optList": optList, "dat": dat, "K": K, "learning_rule": rec_learning_rule,
            #         "hyper": hyper_guess, "weights": weights, "update_w": True, "wMode": None,
            #         "tol": 1e-8, "showOpt": True,
            #        }
            #
            # # Optimization, can also use Nelder-Mead but COBYLA is fastest and pretty reliable
            # res = minimize(evd_lossfun, hyper_to_list(hyper_guess, optList, K), args=args, method='COBYLA')
            # print("Evidence:", -res.fun, "  ", optList, ": ", res.x)
            #
            # opt_hyper = update_hyper(res.x, optList, hyper_guess, K)
            # wMode_sim2, Hess, logEvd, other = getMAP(dat, opt_hyper, weights, W0=None,
            #                                          learning_rule=rec_learning_rule, showOpt=0, tol=1e-12)
            # # Recover error bars for weights
            # W_std = getCredibleInterval(Hess, K)
            # wMode_sim2 = wMode_sim2.reshape((K, -1), order="C")
            #
            # rec_dat = {"args": args, 'res': res, 'opt_hyper': opt_hyper, "W_std": W_std, "wMode": wMode_sim2}
            #
            # fig = plt.figure(figsize=(3.25,1.25))
            # sim_colors = [colors['bias'], colors['h']]
            #
            # for i, c in enumerate(sim_colors):
            #     plt.plot(gen_dat['W'][i], c=c, lw=0.5, zorder=2*i)
            #     plt.plot(rec_dat['wMode'][i], c=c, lw=1, linestyle='--', alpha=0.5, zorder=2*i+1)
            #     plt.fill_between(np.arange(len(rec_dat['wMode'][i])),
            #                      rec_dat['wMode'][i] - 2 * rec_dat['W_std'][i],
            #                      rec_dat['wMode'][i] + 2 * rec_dat['W_std'][i],
            #                      facecolor=c, alpha=0.2, zorder=2*i+1)
            #
            # plt.axhline(0, color="black", linestyle="--", lw=0.5, alpha=0.5, zorder=0)
            #
            # plt.gca().set_xticklabels([])
            # plt.yticks(np.arange(-6,7,2))
            # plt.xlim(0,1000); plt.ylim(-6.2,6.2)
            # # plt.xlabel("Trials"); plt.ylabel("Weights")
            #
            # plt.gca().spines['right'].set_visible(False)
            # plt.gca().spines['top'].set_visible(False)
            #
            # plt.subplots_adjust(0,0,1,1)
            # plt.show()

            # fit to actual data
            #weights = {"bias": 1, "s_a": 1, "s_b": 1}
            weights = {"bias": 1, "cBoth": 1,'stick':1}
            K = np.sum([i for i in weights.values()])

            # save estimated hyper parameters
            if idx==0:
                est_hyper = {
                    'alpha': np.zeros((K,len(ABList))),
                    'sigma':np.zeros((K,len(ABList))),
                    'weight': ['bias', 'cBoth', 'stick'],
                }
            #%%
            # Compute
            hyper_guess = {
                'alpha': [2**-6] * K,
                'sigma': [2**-4] * K,
                'sigInit': [2**4] * K,
                'sigDay': None,
            }

            # Optimizing for both sigma and alpha simultaneously
            optList = ['sigma', 'alpha']

            # List of extra arguments used by evd_lossfun in optimization of evidence
            args = {"optList": optList, "dat": dat, "K": K, "learning_rule": REINFORCE,
                    "hyper": hyper_guess, "weights": weights, "update_w": True, "wMode": None,
                    "tol": 1e-6, "showOpt": True}

            # Optimization, can also use Nelder-Mead but COBYLA is fastest and pretty reliable
            res = minimize(evd_lossfun, hyper_to_list(hyper_guess, optList, K), args=args, method='COBYLA')
            print("Evidence:", -res.fun, "  ", optList, ": ", res.x)

            opt_hyper = update_hyper(res.x, optList, hyper_guess, K)
            wMode, Hess, logEvd, other = getMAP(dat, opt_hyper, weights, W0=None,
                                                learning_rule=rec_learning_rule, showOpt=0, tol=1e-12)
            wMode = wMode.reshape((K, -1), order="C")

            est_hyper['alpha'][:,idx] = opt_hyper['alpha']
            est_hyper['sigma'][:,idx] = opt_hyper['sigma']
            AIC = 2 * K - 2 * logEvd
            BIC = K * np.log(len(dat['y'])) - 2 * logEvd

            AICList[idx,0] = AIC
            BICList[idx,0] = BIC
            negLLList[idx,0] = -logEvd



            # Recover error bars for weights
            #W_std = getCredibleInterval(Hess, K)
            #wMode = wMode.reshape((K, -1), order="C")

            # save

            rec_dat = {"args": args, 'res': res, 'opt_hyper': opt_hyper, "wMode": wMode}

            # recover hyper parameters and std
            hess_args = rec_dat['args'].copy()
            hess_args["wMode"] = rec_dat['wMode'].flatten()
            hess_args["learning_rule"] = hess_args["learning_rule"]
            #H, g = compHess_nolog(evd_lossfun, rec_dat['res'].x, 5e-2, {"keywords": hess_args})
            #hyp_std = np.sqrt(np.diag(np.linalg.inv(H)))

            #rec_dat['hyp_std'] = hyp_std

            # save the file in json
            savedatafolder = os.path.join(dataIndex['BehPath'][dataIndex['Animal']==subject].iloc[0],'latent')
            if not os.path.exists(savedatafolder):
                os.makedirs(savedatafolder)
            savedatapath = os.path.join(dataIndex['BehPath'][dataIndex['Animal']==subject].iloc[0],'latent','psy_fit_'+sessionLabel+'.json')
            with open(savedatapath, 'w') as f:
                json.dump(
                    rec_dat,
                    f,
                    indent=4,
                    default= lambda o: o.tolist() if isinstance(o, np.ndarray)
                                  else int(o) if isinstance(o, np.integer)
                                  else float(o) if isinstance(o, np.floating)
                                  else str(o)
                )

        # # save AIC BIC and log likelihood
        daveFit = {'Animal': subject_inuse, 'AIC': np.squeeze(AICList),
                   'BIC': np.squeeze(BICList), 'negLL': np.squeeze(negLLList)}
        savedatafile = os.path.join(root_dir, strains[0], 'Summary', 'Results', 'PsyFit'+sessionLabel+'.csv')
        pd.DataFrame(daveFit).to_csv(savedatafile)

#%% plot alphas and sigmas by genotype
from scipy.stats import mannwhitneyu
import seaborn as sns
alpha = est_hyper['alpha']
sigma = est_hyper['sigma']
weights = est_hyper['weight']
genotypes = []
for sub in subject_inuse:
    genotypes.append(dataIndex['Genotype'][dataIndex['Animal']==sub].iloc[0])

WTMask = np.array(genotypes) == 'WT'
HetMask = np.array(genotypes) == 'KO'
# Build a tidy DataFrame for plotting
records = []
for i, w in enumerate(weights):
    for subj in range(alpha.shape[1]):
        genotype = 'WT' if WTMask[subj] else 'Het' if HetMask[subj] else 'Other'
        records.append({'param': 'alpha', 'weight': w, 'value': alpha[i, subj], 'genotype': genotype})
        records.append({'param': 'sigma', 'weight': w, 'value': sigma[i, subj], 'genotype': genotype})

df = pd.DataFrame.from_records(records)
plt.figure(figsize=(12,6))
g = sns.catplot(
    data=df, x="weight", y="value", hue="genotype",
    col="param", kind="box", sharey=False,
    height=5, aspect=1
)
for ax, param in zip(g.axes.flat, ['alpha', 'sigma']):
    ymax = df[df['param'] == param]['value'].max()
    for i, w in enumerate(weights):
        subdf = df[(df['param'] == param) & (df['weight'] == w)]
        wt_vals = subdf[subdf['genotype'] == 'WT']['value'].dropna()
        het_vals = subdf[subdf['genotype'] == 'Het']['value'].dropna()

        if len(wt_vals) > 0 and len(het_vals) > 0:
            stat, pval = mannwhitneyu(wt_vals, het_vals, alternative='two-sided')
            # place p-value above the boxes
            ax.text(
                i, ymax * 1.05,
                f"p={pval:.3g}",
                ha='center', va='bottom', fontsize=9, color='black'
            )
    ax.set_title(param.upper())

g.fig.suptitle("Alpha and Sigma by genotype (with Mann–Whitney p-values)", y=1.05)
plt.show()





#%% REINFROCE-base learning rule
# AICList_base = np.zeros(len(ABList))
# BICList_base = np.zeros(len(ABList))
# negLLList_base = np.zeros(len(ABList))
# for idx,ff in enumerate(ABList):
#
#     data = pd.read_csv(ff)
#     #data = data1
#     subject = subjects[idx]
#
# # remove nans
#     data = data.dropna(subset=["actions"])
#     s_a = data['schedule'].copy(); s_b = data['schedule'].copy()
#     s_a[s_a==2] = 0
#     s_b[s_b==1] = 0
#     s_b[s_b==2] = -1
#     correct = data['reward'].copy()
#     correct[np.isnan(correct)] = 0
#     correct[correct>0] = 1
#     answer = (data['schedule'].copy()-1.5)*2
#     y = data['actions'].copy()
#     dayLength = np.array([1])
#
#     # determine inputs
#     t = np.array(data["trial"])
#     prior = ((t[1:] - t[:-1]) == 1).astype(int)
#     prior = np.hstack(([0], prior))
#
#     # Calculate previous average tone value
#     s_avg = (s_a + s_b)
#
#
#     # Calculate previous correct answer
#     h = correct[:-1].astype(int)   # map from (0,1) to (-1,1)
#     h = np.hstack(([0], h))
#     h = h * prior  # for trials without a valid previous trial, set to 0
#
#     # Calculate previous choice
#     c = y[:-1].astype(int)   # map from (0,1) to (-1,1)
#     c = np.hstack(([0], c))
#     c = c * prior  # for trials without a valid previous trial, set to 0
#     stick1 = np.insert(np.array(correct[0:-1])*np.array((y[0:-1]-0.5)*2),0,0)
#
#     stick2 = np.insert(np.array(correct[0:-2])*np.array((y[0:-2]-0.5)*2), 0, [0,0])
#     stick3 = np.insert(np.array(correct[0:-3])*np.array((y[0:-3]-0.5)*2), 0, [0,0,0])
#     inputs = dict(s_a = np.array(s_a)[:, None],
#                   s_b = np.array(s_b)[:, None],
#                   cBoth = np.array(s_avg)[:, None],
#                   stick1 = stick1[:,None],
#                   stick2 = stick2[:,None],
#                   stick3 = stick3[:,None],
#                   h = np.array(h)[:, None],
#                   c = np.array(c)[:, None])
#
#
#     dat = dict(
#         subject = subject,
#         inputs = inputs,
#         s_a = np.array(s_a),
#         s_b = np.array(s_b),
#         correct = np.array(correct),
#         answer = np.array(answer),
#         y =np.array(y),
#         dayLength=dayLength,
#     )
#
#     seed = 42
#     #mouse_name = "CSHL_003"
#     np.random.seed(seed)
#
#     orig_dat = dat
#     #trim_dat = psy.trim(orig_dat, END=N)
#     trim_dat = dat.copy()
#
#
#     rec_learning_rule = REINFORCE
#
#     # fit to actual data
#     #weights = {"bias": 1, "s_a": 1, "s_b": 1}
#     weights = {"bias": 1, "cBoth": 1,'stick1':1, 'stick2':1, 'stick3':1}
#     K = np.sum([i for i in weights.values()])
#
#
#     #%%
#     # Compute
#     hyper_guess = {
#         'alpha': [2**-6] * K,
#         'sigma': [2**-4] * K,
#         'sigInit': [2**4] * K,
#         'sigDay': None,
#     }
#
#     # Optimizing for both sigma and alpha simultaneously
#     optList = ['sigma', 'alpha']
#
#     # List of extra arguments used by evd_lossfun in optimization of evidence
#     args = {"optList": optList, "dat": dat, "K": K, "learning_rule": REINFORCE,
#             "hyper": hyper_guess, "weights": weights, "update_w": True, "wMode": None,
#             "tol": 1e-6, "showOpt": True}
#
#     # Optimization, can also use Nelder-Mead but COBYLA is fastest and pretty reliable
#     res = minimize(evd_lossfun, hyper_to_list(hyper_guess, optList, K), args=args, method='COBYLA')
#     print("Evidence:", -res.fun, "  ", optList, ": ", res.x)
#
#     opt_hyper = update_hyper(res.x, optList, hyper_guess, K)
#     wMode, Hess, logEvd, other = getMAP(dat, opt_hyper, weights, W0=None,
#                                         learning_rule=rec_learning_rule, showOpt=0, tol=1e-12)
#
#     AIC = 2 * K - 2 * logEvd
#     BIC = K * np.log(len(dat['y'])) - 2 * logEvd
#
#     AICList_base[idx] = AIC
#     BICList_base[idx] = BIC
#     negLLList_base[idx] = -logEvd
#
#     # Recover error bars for weights
#     W_std = getCredibleInterval(Hess, K)
#     wMode = wMode.reshape((K, -1), order="C")
#
#     rec_dat = {"args": args, 'res': res, 'opt_hyper': opt_hyper, "W_std": W_std, "wMode": wMode}
#
#



# %%
sim_colors = [colors['bias'], colors['cBoth'], colors['c']]

plt.figure(figsize=(3.25, 1.25))
for i, c in enumerate(sim_colors):
    plt.plot(wMode[i], c=c, lw=1, linestyle='-', alpha=0.85, zorder=2 * i + 1)
    plt.fill_between(np.arange(len(wMode[i])), wMode[i] - 2 * W_std[i], wMode[i] + 2 * W_std[i],
                     facecolor=c, alpha=0.2, zorder=2 * i)

plt.axhline(0, color="black", linestyle="--", lw=0.5, alpha=0.5, zorder=0)
plt.xticks(1000 * np.arange(0, 11))
plt.yticks(np.arange(-2, 3, 1))
plt.xlim(0, 1000)
plt.ylim(-2.25, 2.25)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.show()

#%% plot recoverd hyperparameters
plt.figure(figsize=(1.5, 1.25))
sim_colors = [colors['bias'], colors['h']]
num_std = 1.00
K = gen_dat['K']

true_sigma = gen_dat['true_sigma']
avg_sigma = np.log2(rec_dat['opt_hyper']['sigma'])
err_sigma = rec_dat['hyp_std'][:K]
for i, c in enumerate(sim_colors):
    yerr = [[-(np.log2(2 ** avg_sigma[i] - err_sigma[i] * num_std) - avg_sigma[i])],
            [np.log2(2 ** avg_sigma[i] + err_sigma[i] * num_std) - avg_sigma[i]]]
    if np.isnan(yerr[0][0]):
        yerr[0][0] = 100
    plt.plot([i - 0.3, i + 0.3], [np.log2(true_sigma[i])] * 2, color=c, linestyle="-", lw=1.2, zorder=1, alpha=0.8)
    plt.errorbar([i], avg_sigma[i], yerr=yerr, c=c, lw=1, marker='o', markersize=4)

true_alpha = gen_dat['true_alpha']
avg_alpha = np.log2(rec_dat['opt_hyper']['alpha'])
err_alpha = rec_dat['hyp_std'][K:2 * K]
for i, c in enumerate(sim_colors):
    yerr = [[-(np.log2(2 ** avg_alpha[i] - err_alpha[i] * num_std) - avg_alpha[i])],
            [np.log2(2 ** avg_alpha[i] + err_alpha[i] * num_std) - avg_alpha[i]]]
    if np.isnan(yerr[0][0]):
        yerr[0][0] = 100
    plt.plot([i + 1.7, i + 2.3], [np.log2(true_alpha[i])] * 2, color=c, linestyle="-", lw=1.2, zorder=1, alpha=0.8)
    plt.errorbar([i + 2], avg_alpha[i], yerr=yerr, c=c, lw=1, marker='s', markersize=4)

plt.ylim(-8.75, -3.5)
plt.yticks(np.arange(-8, -3))
plt.xticks([0, 1, 2, 3])
plt.gca().set_xticklabels([r"$\sigma_1$", r"$\sigma_2$",
                           r"$\alpha_1$", r"$\alpha_2$"])

# %%
# %%
v_rw = other['pT']['learning_terms']['v'].T
E_rw = other['pT']['learning_terms']["E_flat"].reshape((K, -1), order="C")

noise = np.cumsum(E_rw, axis=1)
learning = np.cumsum(v_rw, axis=1)

plt.figure(figsize=(3.25, 1.25))
for i, c in enumerate(sim_colors):
    plt.plot(learning[i], c=c, lw=0.75, linestyle='-', alpha=0.85, zorder=2 * i)
    plt.fill_between(np.arange(len(noise[i])), learning[i], learning[i] + noise[i],
                     facecolor=c, alpha=0.2, zorder=2 * i)
    plt.plot(learning[i] + noise[i], c=c, lw=0.25, linestyle='-', alpha=0.5, zorder=2 * i + 1)

plt.axhline(0, color="black", linestyle="--", lw=0.5, alpha=0.5, zorder=0)
plt.xticks(1000 * np.arange(0, 11))
plt.yticks(np.arange(-2, 3, 1))
plt.xlim(0, 1000)
plt.ylim(-2.25, 2.25)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.show()

#%%
hess_args = rec_dat['args'].copy()
hess_args["wMode"] = rec_dat['wMode'].flatten()
hess_args["learning_rule"] = hess_args["learning_rule"]
H, g = compHess_nolog(evd_lossfun, rec_dat['res'].x, 5e-2, {"keywords": hess_args})
hyp_std = np.sqrt(np.diag(np.linalg.inv(H)))

rec_dat['hyp_std'] = hyp_std

# Plot recovered hyperparameters
